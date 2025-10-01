---
output: rmarkdown::github_document
---

<!-- README.md is generated from README.Rmd. Please edit that file -->



# bigPCAcpp

# Principal component analysis for bigmemory matrices
## Frédéric Bertrand

<!-- badges: start -->
[![R-CMD-check](https://github.com/fbertran/bigPCAcpp/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/fbertran/bigPCAcpp/actions/workflows/R-CMD-check.yaml)
<!-- badges: end -->

The **bigPCAcpp** package provides high performance principal component
analysis (PCA) routines for `bigmemory::big.matrix` objects. The package
streams data through BLAS and LAPACK kernels so that very large matrices
can be analysed without first copying them into native R matrices. It
also includes helpers to export PCA diagnostics such as scores,
loadings, correlations, and contributions into file-backed
`big.matrix` allocations suitable for downstream analysis pipelines.

`bigPCAcpp` exposes two complementary backends:

* **bigmemory** – operates directly on in-memory or file-backed
  `big.matrix` objects using the bigmemory C++ API.
* **bigalgebra** – writes results into `bigalgebra` matrices so that
  large studies can keep intermediate artifacts out-of-memory while
  benefiting from native BLAS/LAPACK performance.

These workflows make it possible to analyse data sets that exceed the
available RAM while keeping numerical stability through double-precision
accumulation and LAPACK eigen decompositions.

## Installation

You can install the development version of **bigPCAcpp** from GitHub with:


``` r
# install.packages("devtools")
devtools::install_github("fbertran/bigPCAcpp")
```

If you prefer a local source install, clone the repository and run:


``` r
R CMD build bigPCAcpp
R CMD INSTALL bigPCAcpp_0.0.1.tar.gz
```

For workflows that rely on the streaming big.matrix helpers, ensure that the
`bigalgebra` package is installed:


``` r
install.packages("bigalgebra")
```

## Options

The package defines several options to control numerical tolerances and
workspace allocation. They are prefixed with `bigPCAcpp.` and include:

Option | Default value | Description
--- | --- | ---
`bigPCAcpp.block_size` | `1000L` | Number of rows processed in each block when streaming scores through BLAS.
`bigPCAcpp.center_scale_epsilon` | `1e-8` | Lower bound applied when rescaling columns to avoid division instabilities.
`bigPCAcpp.progress` | `FALSE` | Emit progress updates when computing PCA on long-running jobs.

All options can be changed with `options()` at runtime. For example,
`options(bigPCAcpp.block_size = 5000L)` increases the streaming block size.

## Examples

The examples below demonstrate the bigmemory workflow and compare the
results with base R's `prcomp()` implementation.


``` r
library(bigmemory)
library(bigPCAcpp)

# Allocate a 1,000 x 25 big.matrix with simulated values
n <- 1000
p <- 25
bm <- bigmemory::big.matrix(n, p, type = "double")
bm[,] <- matrix(rnorm(n * p), nrow = n)

# Run PCA and extract eigenvalues and rotation
res <- pca_bigmatrix(bm, center = TRUE, scale. = TRUE)
res$importance
res$rotation[1:5, 1:3]

# Generate PCA scores in bigmemory-backed storage
scores <- bigmemory::filebacked.big.matrix(n, ncol = 3, type = "double")
pca_scores_stream_bigmatrix(
  bm,
  scores,
  res$rotation[, 1:3],
  center = res$center,
  scale = res$scale
)
scores[1:5, ]

# Compare with prcomp()
pr <- prcomp(bm[], center = TRUE, scale. = TRUE)
all.equal(unclass(pr$rotation)[, 1:3], res$rotation[, 1:3], tolerance = 1e-6)
```

`pca_bigmatrix()` can also focus on a subset of leading components while
streaming the results into file-backed matrices. The following snippet stores the
first four principal components and keeps a running summary of their scores.


``` r
library(bigmemory)
library(bigPCAcpp)

set.seed(2025)
bm <- bigmemory::filebacked.big.matrix(nrow = 1500, ncol = 40, type = "double")
bm[,] <- matrix(rnorm(1500 * 40), nrow = 1500)

# Request only the first four components
top_pca <- pca_bigmatrix(bm, center = TRUE, scale. = TRUE, ncomp = 4)
top_pca$sdev

# Stream the corresponding scores into a file-backed allocation
scores_fb <- bigmemory::filebacked.big.matrix(nrow = nrow(bm), ncol = 4, type = "double")
pca_scores_stream_bigmatrix(
  bm,
  scores_fb,
  top_pca$rotation[, 1:4],
  center = top_pca$center,
  scale = top_pca$scale
)

# Inspect a lightweight summary without loading the entire matrix
colMeans(scores_fb[, 1:2])
apply(scores_fb[, 1:2], 2, sd)
```

To stream the diagnostics into `bigmemory`-backed matrices, use the
corresponding helper functions:


``` r
library(bigmemory)
library(bigalgebra)
library(bigPCAcpp)

n <- 1000
p <- 25
bm <- bigmemory::filebacked.big.matrix(n, p, type = "double")
bm[,] <- matrix(rnorm(n * p), nrow = n)

rotation <- bigmemory::big.matrix(nrow = p, ncol = p)
loadings <- bigmemory::big.matrix(nrow = p, ncol = p)
correlations <- bigmemory::big.matrix(nrow = p, ncol = p)
contrib <- bigmemory::big.matrix(nrow = p, ncol = p)

pca_stream <- pca_stream_bigmatrix(bm, xpRotation = rotation, center = TRUE, scale = FALSE)
pca_variable_loadings_stream_bigmatrix(rotation, pca_stream$sdev, loadings)
pca_variable_correlations_stream_bigmatrix(rotation, pca_stream$sdev, correlations)
pca_variable_contributions_stream_bigmatrix(loadings, contrib)
```

### Robust PCA and singular value decompositions

Robust workflows dampen the influence of outliers while retaining the
familiar PCA interface. The `pca_robust()` helper centres variables by the
median, optionally scales by the MAD, and relies on an iteratively
reweighted SVD to derive principal components. The same robust solver is
exposed directly via `svd_robust()` for use in custom pipelines, and the
streaming-friendly `svd_bigmatrix()` wrapper computes classical SVDs on
`big.matrix` objects without materialising dense copies in memory.


``` r
library(bigmemory)
library(bigPCAcpp)

set.seed(42)
mat <- matrix(rnorm(200), nrow = 40, ncol = 5)
mat[1, 1] <- 15  # introduce an outlier

# Robust PCA keeps the outlier from dominating the rotation
robust <- pca_robust(mat, ncomp = 3)
robust$explained_variance

# Classical PCA on the same data highlights the impact of the outlier
bm_small <- bigmemory::big.matrix(nrow = nrow(mat), ncol = ncol(mat), type = "double")
bm_small[,] <- mat
classical <- pca_bigmatrix(bm_small, center = TRUE, scale. = TRUE, ncomp = 3)
cbind(classical = classical$rotation[1:5, 1], robust = robust$rotation[1:5, 1])

# Classical SVD on a file-backed big.matrix
bm <- bigmemory::filebacked.big.matrix(200, 10, type = "double")
bm[,] <- matrix(rnorm(2000), nrow = 200)
svd_stream <- svd_bigmatrix(bm, nu = 3, nv = 3)
svd_stream$d

# Direct access to the robust SVD routine
svd_out <- svd_robust(mat, ncomp = 3)
svd_out$d
svd_out$weights[1:6]
```

Robust decompositions down-weight the contaminated observations while the
classical stream demonstrates how to fetch singular vectors without materialising
the dense matrix. The robust solver also exposes per-row weights that can be
reused to flag problematic observations for further inspection.

### Plotting diagnostics

`bigPCAcpp` bundles plot helpers that operate on both dense matrices and
`big.matrix` backends. The snippets below illustrate how to call each
function using results from `pca_bigmatrix()`. For instance, the
`pca_plot_scores()` helper samples observations and draws a scatter plot of
their scores on a chosen pair of components, which is particularly useful when
you need to visually assess potential clusters without loading the full data
set into memory.


``` r
library(bigmemory)
library(bigPCAcpp)

set.seed(123)
bm <- bigmemory::big.matrix(500, 6, type = "double")
bm[,] <- matrix(rnorm(500 * 6), nrow = 500)
res <- pca_bigmatrix(bm, center = TRUE, scale. = TRUE)

# Scree plot of explained variance
pca_plot_scree(res)

# Scatter plot of sampled scores on PCs 1 and 2
pca_plot_scores(
  bm,
  res$rotation,
  center = res$center,
  scale = res$scale,
  components = c(1L, 2L),
  max_points = 2000L,
  seed = 2024
)

# Contribution bar plot for the leading component
loadings <- pca_variable_loadings(res$rotation, res$sdev)
contrib <- pca_variable_contributions(loadings)
pca_plot_contributions(contrib, component = 1L, top_n = 10L)

# Correlation circle for the first two components
correlations <- pca_variable_correlations(res$rotation, res$sdev, res$column_sd)
pca_plot_correlation_circle(correlations, components = c(1L, 2L))

# Biplot combining scores and loadings
scores <- res$scores
if (is.null(scores)) {
  scores <- pca_scores_bigmatrix(bm, res$rotation, center = res$center, scale = res$scale)
}
pca_plot_biplot(scores, loadings, components = c(1L, 2L))
```

## Citation

If you use **bigPCAcpp** in academic work, please cite:

Bertrand F. (2025). *bigPCAcpp: Principal Component Analysis for bigmemory Matrices*.

## Maintainer

Maintainer: Frédéric Bertrand <frederic.bertrand@lecnam.net>

For questions, bug reports, or contributions, please open an issue on
[GitHub](https://github.com/fbertran/bigPCAcpp).

