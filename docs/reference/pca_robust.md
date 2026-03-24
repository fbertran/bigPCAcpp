# Robust principal component analysis

Compute principal component analysis (PCA) using robust measures of
location and scale so that extreme observations have a reduced influence
on the resulting components. The implementation centres each variable by
its median and, when requested, scales by the median absolute deviation
(MAD) before performing an iteratively reweighted singular value
decomposition that down-weights observations with unusually large
reconstruction errors.

## Usage

``` r
pca_robust(x, center = TRUE, scale = FALSE, ncomp = NULL)
```

## Arguments

- x:

  A numeric matrix, data frame, or an object coercible to a numeric
  matrix. Missing values are not supported.

- center:

  Logical; should variables be centred by their median before applying
  PCA?

- scale:

  Logical; when `TRUE`, variables are scaled by the MAD after centring.
  Scaling requires `center = TRUE`.

- ncomp:

  Number of components to retain. Use `NULL` or a non-positive value to
  keep all components returned by the decomposition.

## Value

A [`bigpca`](https://fbertran.github.io/bigPCAcpp/reference/bigpca.md)
object mirroring the structure of
[`pca_bigmatrix()`](https://fbertran.github.io/bigPCAcpp/reference/pca_bigmatrix.md)
with robust estimates of location, scale, and variance metrics.

## Examples

``` r
set.seed(42)
x <- matrix(rnorm(50), nrow = 10)
x[1, 1] <- 25  # outlier
robust <- pca_robust(x, ncomp = 2)
robust$sdev
#> [1] 2.801017 1.603985
```
