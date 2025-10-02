---
output: rmarkdown::github_document
---

<!-- README.md is generated from README.Rmd. Please edit that file -->



# bigPCAcpp <img src="man/figures/BigPCAcpp_hex_dark.svg" align="right" width="200"/>

# Principal component analysis for bigmemory matrices
## Frédéric Bertrand

<!-- badges: start -->
[![R-CMD-check](https://github.com/fbertran/bigPCAcpp/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/fbertran/bigPCAcpp/actions/workflows/R-CMD-check.yaml)
[![R-hub](https://github.com/fbertran/bigPCAcpp/actions/workflows/rhub.yaml/badge.svg)](https://github.com/fbertran/bigPCAcpp/actions/workflows/rhub.yaml)
<!-- badges: end -->

The **bigPCAcpp** package provides high performance principal component
analysis (PCA) routines specialised for `bigmemory::big.matrix` objects.
It keeps data in bigmemory allocations from ingestion through
eigendecomposition so that very large matrices can be analysed without
copying them into base R matrices. In addition to the PCA core, the
package offers streaming helpers that write scores, loadings,
correlations, and contributions back into file-backed `big.matrix`
targets for integration with downstream pipelines.

Beyond classical PCA, the package ships with scalable SVD tools that can
process file-backed matrices block by block, and it includes robust PCA
and robust SVD routines that temper the influence of outliers while
remaining compatible with bigmemory workflows. For exploratory work on
large batches, a scalable PCA interface lets users extract leading
components without reading the full matrix into memory.

These workflows make it possible to analyse data sets that exceed the
available RAM while keeping numerical stability through double-precision
accumulation and LAPACK eigen decompositions. Current features include

* centring and scaling directly on `big.matrix` inputs,
* incremental score generation that avoids bringing data into memory,
* helpers to persist PCA diagnostics in file-backed storage,
* SVD utilities for dense and robust decompositions backed by
  `bigmemory`, and

## Installation

You can install the development version of **bigPCAcpp** from GitHub with:


``` r
# install.packages("devtools")
devtools::install_github("fbertran/bigPCAcpp")
```

If you prefer a local source install, clone the repository and run:


``` r
R CMD build bigPCAcpp
R CMD INSTALL bigPCAcpp_0.9.0.tar.gz
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
res <- pca_bigmatrix(bm, center = TRUE, scale = TRUE)
res$eigenvalues
#>  [1] 1.2701077 1.2537526 1.2151939 1.2024274 1.1802276 1.1354639 1.1048690
#>  [8] 1.0902075 1.0630117 1.0403948 1.0224523 1.0122884 0.9931732 0.9852838
#> [15] 0.9550143 0.9414302 0.9204860 0.9008016 0.8898188 0.8773914 0.8437966
#> [22] 0.7979435 0.7821215 0.7628152 0.7595271
res$importance
#> NULL
res$rotation[1:5, 1:3]
#>             [,1]        [,2]        [,3]
#> [1,]  0.14398208  0.22870330  0.12242692
#> [2,] -0.07237026 -0.21951399  0.19761877
#> [3,] -0.32741500 -0.34723923 -0.02282739
#> [4,]  0.23671229  0.14764528  0.05672112
#> [5,]  0.14343191  0.03490776  0.10010376

# Generate PCA scores in bigmemory storage
scores <- bigmemory::big.matrix(
  nrow = n,
  ncol = 3,
  type = "double"
)

(pca_scores_bigmatrix(
  bm,
  res$rotation,
  center = res$center,
  scale = res$scale
))[1:6,1:6]
#>             [,1]       [,2]       [,3]       [,4]        [,5]         [,6]
#> [1,]  0.65907229 -0.1974698 -3.5714713 -0.7150943  0.01412313  0.519975663
#> [2,] -0.20615917  2.6536525 -2.1246929 -0.2443652 -1.22009504  0.338157676
#> [3,]  1.97067655 -0.1459897  0.2706729 -2.5778456  1.03548532  1.343923490
#> [4,] -0.06521853 -0.3507049 -0.7159949  0.2779695  1.08616408 -1.106634062
#> [5,]  0.34681884  0.6974558 -2.4466572 -0.4696698  1.95794498 -0.008496806
#> [6,]  0.62390402 -0.9625778 -0.9364008  1.5785386  0.47807633  0.074748751

# Compare sum of absolute values with prcomp()
pr <- prcomp(bm[], center = TRUE, scale = TRUE)
sum(abs(abs(pr$rotation[, 1:3])-abs(res$rotation[, 1:3])))<10^(-6)
#> [1] TRUE
```

`pca_bigmatrix()` can also focus on a subset of leading components while
streaming the results into file-backed matrices. The following snippet stores the
first four principal components and keeps a running summary of their scores.


``` r
library(bigmemory)
library(bigPCAcpp)

set.seed(2025)
bm <- bigmemory::big.matrix(nrow = 1500, ncol = 40, type = "double")
bm[,] <- matrix(rnorm(1500 * 40), nrow = 1500)

# Request only the first four components
top_pca <- pca_bigmatrix(bm, center = TRUE, scale = TRUE, ncomp = 4)
top_pca$sdev
#> [1] 1.141546 1.124998 1.119607 1.109924

# Stream the corresponding scores into a file-backed allocation
path <- tempfile(fileext = ".bin")
desc <- paste0(path, ".desc")

scores_fb <- bigmemory::filebacked.big.matrix(nrow = nrow(bm), ncol = 4,
             type = "double", backingfile = basename(path), backingpath =
             dirname(path), descriptorfile = basename(desc)
)
pca_scores_stream_bigmatrix(
  bm,
  scores_fb,
  top_pca$rotation[, 1:4],
  center = top_pca$center,
  scale = top_pca$scale
)
#> <pointer: 0x1596ceca0>

# Inspect a lightweight summary without loading the entire matrix
colMeans(scores_fb[, 1:2])
#> [1] 3.944992e-17 3.064216e-17
apply(scores_fb[, 1:2], 2, sd)
#> [1] 1.141546 1.124998
```

To stream the diagnostics into `bigmemory`-backed matrices, use the
corresponding helper functions:


``` r
library(bigmemory)
library(bigPCAcpp)

n <- 1000
p <- 25
bm <- bigmemory::big.matrix(n, p, type = "double")
bm[,] <- matrix(rnorm(n * p), nrow = n)

rotation <- bigmemory::big.matrix(nrow = p, ncol = p)
loadings <- bigmemory::big.matrix(nrow = p, ncol = p)
correlations <- bigmemory::big.matrix(nrow = p, ncol = p)
contrib <- bigmemory::big.matrix(nrow = p, ncol = p)

pca_stream <- pca_stream_bigmatrix(bm, xpRotation = rotation, 
                                   center = TRUE, scale = FALSE)
pca_variable_loadings_stream_bigmatrix(rotation, pca_stream$sdev,
                                       loadings)
#> <pointer: 0x1596df4e0>
pca_variable_correlations_stream_bigmatrix(rotation, pca_stream$sdev,
                           pca_stream$column_sd, correlations)
#> Error in pca_variable_correlations_stream_bigmatrix(rotation, pca_stream$sdev, : argument "xpDest" is missing, with no default
pca_variable_contributions_stream_bigmatrix(loadings, contrib)
#> <pointer: 0x159655010>
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
mat_scaled <- scale(mat, center = TRUE, scale=TRUE)

# Classical PCA on the same data highlights the impact of the outlier
bm_small <- bigmemory::big.matrix(nrow = nrow(mat_scaled), ncol = ncol(mat_scaled), type = "double")
bm_small[,] <- mat_scaled
classical <- pca_bigmatrix(bm_small, center = FALSE, scale = FALSE, ncomp = 3)
classical$explained_variance
#> [1] 0.2940708 0.2332728 0.2031007

scores_classical <- pca_scores_bigmatrix(xpMat = bm_small, rotation = classical$rotation, center = classical$center, classical$scale)
scores_classical[1,]
#> [1] -4.752614 -1.534966  1.578737

pca_plot_contributions(pca_individual_contributions(scores_classical, classical$sdev))
```

<div class="figure">
<img src="man/figures/README-robustsvdexample-1.png" alt="plot of chunk robustsvdexample" width="100%" />
<p class="caption">plot of chunk robustsvdexample</p>
</div>

``` r

# Robust PCA keeps the outlier from dominating the rotation
robust <- pca_robust(mat_scaled, center = FALSE, scale = FALSE, ncomp = 3)
robust$explained_variance
#> [1] 0.3633363 0.3509611 0.2857026

robust$scores[1,]
#> [1] 1.025663 1.948710 2.095546

pca_plot_contributions(pca_individual_contributions(robust$scores, robust$sdev))
```

<div class="figure">
<img src="man/figures/README-robustsvdexample-2.png" alt="plot of chunk robustsvdexample" width="100%" />
<p class="caption">plot of chunk robustsvdexample</p>
</div>

``` r

cbind(classical = classical$rotation[1:5, 1], robust = robust$rotation[1:5, 1])
#>       classical     robust
#> [1,] -0.5793644 0.01128235
#> [2,] -0.3121420 0.59547597
#> [3,] -0.5716000 0.77399456
#> [4,]  0.4071138 0.18028142
#> [5,]  0.2728298 0.11709868
```


``` r
# Classical SVD on a file-backed big.matrix
path <- tempfile(fileext = ".bin")
desc <- paste0(path, ".desc")

bm <- bigmemory::filebacked.big.matrix(200, 10, type = "double", backingfile = 
      basename(path), backingpath = dirname(path), descriptorfile = basename(desc))
bm[,] <- matrix(rnorm(2000), nrow = 200)
svd_stream <- svd_bigmatrix(bm, nu = 3, nv = 3)
svd_stream$d
#>  [1] 16.66256 15.90085 15.80823 14.84659 13.99062 13.52699 13.06717 12.61343
#>  [9] 12.15871 11.63997

# Direct access to the robust SVD routine
svd_out <- svd_robust(mat, ncomp = 3)
svd_out$d
#> [1] 16.789433  6.178555  5.620833
svd_out$weights[1:6]
#> [1] 1 1 1 1 1 1
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
res <- pca_bigmatrix(bm, center = TRUE, scale = TRUE)

# Scree plot of explained variance
pca_plot_scree(res)
```

<div class="figure">
<img src="man/figures/README-plotexamples-1.png" alt="plot of chunk plotexamples" width="100%" />
<p class="caption">plot of chunk plotexamples</p>
</div>

``` r

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
```

<div class="figure">
<img src="man/figures/README-plotexamples-2.png" alt="plot of chunk plotexamples" width="100%" />
<p class="caption">plot of chunk plotexamples</p>
</div>

``` r

# Contribution bar plot for the leading component
loadings <- pca_variable_loadings(res$rotation, res$sdev)
contrib <- pca_variable_contributions(loadings)
pca_plot_contributions(contrib, component = 1L, top_n = 10L)
```

<div class="figure">
<img src="man/figures/README-plotexamples-3.png" alt="plot of chunk plotexamples" width="100%" />
<p class="caption">plot of chunk plotexamples</p>
</div>

``` r

# Correlation circle for the first two components
correlations <- pca_variable_correlations(res$rotation, res$sdev, 
res$column_sd, res$scale)
pca_plot_correlation_circle(correlations, components = c(1L, 2L))
```

<div class="figure">
<img src="man/figures/README-plotexamples-4.png" alt="plot of chunk plotexamples" width="100%" />
<p class="caption">plot of chunk plotexamples</p>
</div>

``` r

# Biplot combining scores and loadings
scores <- res$scores
if (is.null(scores)) {
  scores <- pca_scores_bigmatrix(bm, res$rotation, center = res$center, scale = res$scale)
}
pca_plot_biplot(scores, loadings, components = c(1L, 2L))
```

<div class="figure">
<img src="man/figures/README-plotexamples-5.png" alt="plot of chunk plotexamples" width="100%" />
<p class="caption">plot of chunk plotexamples</p>
</div>

## Citation

If you use **bigPCAcpp** in academic work, please cite:

Bertrand F. (2025). *bigPCAcpp: Principal Component Analysis for bigmemory Matrices*.

## Maintainer

Maintainer: Frédéric Bertrand <frederic.bertrand@lecnam.net>

For questions, bug reports, or contributions, please open an issue on
[GitHub](https://github.com/fbertran/bigPCAcpp).

