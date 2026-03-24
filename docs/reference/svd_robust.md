# Robust singular value decomposition (C++ backend)

Compute the iteratively reweighted SVD using the high-performance C++
implementation. The interface mirrors
[`svd_robust_R()`](https://fbertran.github.io/bigPCAcpp/reference/svd_robust_R.md)
while delegating the heavy lifting to compiled code.

## Usage

``` r
svd_robust(
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
of iterations required for convergence (`iterations`).
