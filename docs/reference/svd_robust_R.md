# Iteratively reweighted singular value decomposition

Internal helper used by
[`pca_robust()`](https://fbertran.github.io/bigPCAcpp/reference/pca_robust.md)
to compute a singular value decomposition that is less sensitive to
individual rows with extreme values. The routine alternates between
computing the SVD of a row-weighted matrix and updating the weights via
a Huber-type scheme based on the reconstruction residuals.

## Usage

``` r
svd_robust_R(
  x,
  ncomp,
  max_iter = 25L,
  tol = sqrt(.Machine$double.eps),
  huber_k = 1.345
)
```

## Arguments

- x:

  Numeric matrix for which the decomposition should be computed.

- ncomp:

  Number of leading components to retain.

- max_iter:

  Maximum number of reweighting iterations.

- tol:

  Convergence tolerance applied to successive changes in the row weights
  and singular values.

- huber_k:

  Tuning constant controlling the aggressiveness of the Huber weight
  function. Larger values down-weight fewer observations.

## Value

A list containing the left and right singular vectors (`u` and `v`), the
singular values (`d`), the final row weights (`weights`), and the number
of iterations required for convergence (`iterations`). The structure
mirrors base R's [`base::svd()`](https://rdrr.io/r/base/svd.html) output
with additional metadata.
