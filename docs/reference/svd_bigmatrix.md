# Singular value decomposition for `bigmemory::big.matrix` inputs

Compute the singular value decomposition (SVD) of a
[`bigmemory::big.matrix`](https://rdrr.io/pkg/bigmemory/man/big.matrix.html)
without materialising it as a base R matrix. Blocks of rows are streamed
through BLAS before LAPACK is invoked so that even moderately large
matrices can be decomposed efficiently.

## Usage

``` r
svd_bigmatrix(
  xpMat,
  nu = -1L,
  nv = -1L,
  block_size = 1024L,
  method = c("dgesdd", "dgesvd")
)
```

## Arguments

- xpMat:

  Either a
  [`bigmemory::big.matrix`](https://rdrr.io/pkg/bigmemory/man/big.matrix.html)
  or an external pointer such as `mat@address` that references the
  source `big.matrix`.

- nu:

  Number of left singular vectors to return. Use a negative value to
  request the default of `min(nrow, ncol)` vectors and zero to skip
  returning `u` entirely.

- nv:

  Number of right singular vectors to return. Use a negative value to
  request the default of `min(nrow, ncol)` vectors and zero to skip
  returning `v` entirely.

- block_size:

  Number of rows to process per block when streaming data into BLAS
  kernels. Larger values can improve throughput at the cost of
  additional temporary memory.

- method:

  LAPACK backend used to compute the decomposition. The default uses the
  divide-and-conquer routine `dgesdd` and falls back to `dgesvd` when
  required.

## Value

A list with components `u`, `d`, and `v` analogous to base R's
[`svd()`](https://rdrr.io/r/base/svd.html) output. When `nu` or `nv` are
zero the corresponding matrix has zero columns.

## Examples

``` r
set.seed(42)
mat <- bigmemory::as.big.matrix(matrix(rnorm(20), nrow = 5))
svd_res <- svd_bigmatrix(mat, nu = 2, nv = 2)
svd_res$d
#> [1] 4.319302 3.283675 1.685888 1.091984
```
