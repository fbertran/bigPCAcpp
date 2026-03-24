# Prepare iteratively reweighted singular value decomposition

Internal helper used by
[`pca_robust()`](https://fbertran.github.io/bigPCAcpp/reference/pca_robust.md)
to compute a singular value decomposition that is less sensitive to
individual rows with extreme values. The routine alternates between
computing the SVD of a row-weighted matrix and updating the weights via
a Huber-type scheme based on the reconstruction residuals.

## Usage

``` r
prepare_svd_robust_input(x, ncomp, max_iter, tol, huber_k)
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

A list containing x, n, p, ncomp, max_iter, tol and huber_k.
