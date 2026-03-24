# Streaming big.matrix PCA helpers

Variants of the PCA helpers that stream results directly into
[`bigmemory::big.matrix`](https://rdrr.io/pkg/bigmemory/man/big.matrix.html)
objects, enabling file-backed workflows without materialising dense R
matrices.

## Usage

``` r
pca_spca_stream_bigmatrix(
  xpMat,
  xpRotation = NULL,
  center = TRUE,
  scale = FALSE,
  ncomp = -1L,
  block_size = 2048L,
  max_iter = 50L,
  tol = 1e-04,
  seed = NULL,
  return_scores = FALSE,
  verbose = FALSE
)

pca_scores_stream_bigmatrix(
  xpMat,
  xpDest,
  rotation,
  center,
  scale,
  ncomp = -1L,
  block_size = 1024L
)

pca_variable_loadings_stream_bigmatrix(xpRotation, sdev, xpDest)

pca_variable_correlations_stream_bigmatrix(
  xpRotation,
  sdev,
  column_sd,
  scale = NULL,
  xpDest
)

pca_variable_contributions_stream_bigmatrix(xpLoadings, xpDest)
```

## Arguments

- xpMat:

  Either a
  [`bigmemory::big.matrix`](https://rdrr.io/pkg/bigmemory/man/big.matrix.html)
  or an external pointer such as `mat@address` that references the
  source `big.matrix`.

- xpRotation:

  For `pca_variable_correlations_stream_bigmatrix()`, a
  [`bigmemory::big.matrix`](https://rdrr.io/pkg/bigmemory/man/big.matrix.html)
  or external pointer containing the rotation matrix to stream from.

- center:

  For
  [`pca_scores_bigmatrix()`](https://fbertran.github.io/bigPCAcpp/reference/pca_bigmatrix.md),
  a numeric vector of column means (optional).

- scale:

  Optional numeric vector of scaling factors returned by
  `pca_stream_bigmatrix()` or
  [`pca_bigmatrix()`](https://fbertran.github.io/bigPCAcpp/reference/pca_bigmatrix.md).
  When supplied, correlations are reported on the scaled data without
  dividing by `column_sd`.

- ncomp:

  Number of components to retain. Use a non-positive value to keep all
  components returned by the decomposition.

- block_size:

  Number of rows to process per block when streaming data through BLAS
  kernels. Larger values improve throughput at the cost of additional
  memory.

- max_iter:

  Maximum number of block power iterations.

- tol:

  Convergence tolerance applied to the Frobenius norm of the difference
  between successive subspace projectors.

- seed:

  Optional integer seed used to initialise the random starting basis.

- return_scores:

  Logical; when `TRUE`, principal component scores are computed in a
  final streaming pass over the data.

- verbose:

  Logical; when `TRUE`, diagnostic messages describing the iteration
  progress are emitted.

- xpDest:

  Either a `big.matrix` or external pointer referencing the destination
  `big.matrix` that stores the computed quantity.

- rotation:

  A rotation matrix such as the `rotation` element returned by
  [`pca_bigmatrix()`](https://fbertran.github.io/bigPCAcpp/reference/pca_bigmatrix.md).

- sdev:

  A numeric vector of component standard deviations, typically the
  `sdev` element from
  [`pca_bigmatrix()`](https://fbertran.github.io/bigPCAcpp/reference/pca_bigmatrix.md).

- column_sd:

  A numeric vector of variable standard deviations used to scale the
  correlations when the PCA was performed on unscaled data.

- xpLoadings:

  For `pca_variable_contributions_stream_bigmatrix()`, the loadings
  matrix supplied as a `big.matrix` or external pointer.

## Value

For `pca_stream_bigmatrix()`, the same `bigpca` object as
[`pca_bigmatrix()`](https://fbertran.github.io/bigPCAcpp/reference/pca_bigmatrix.md)
with the addition of a `rotation_stream_bigmatrix` element referencing
the populated `big.matrix` when `xpRotation` is supplied. For
`pca_spca_stream_bigmatrix()`, the same scalable PCA structure as
[`pca_spca()`](https://fbertran.github.io/bigPCAcpp/reference/pca_spca.md)
with the optional pointer populated when provided.

The external pointer supplied in `xpDest`, invisibly.

## Functions

- `pca_scores_stream_bigmatrix()`: Stream PCA scores into a destination
  big.matrix.

- `pca_variable_loadings_stream_bigmatrix()`: Populate big.matrix
  objects with derived variable diagnostics.

- `pca_variable_correlations_stream_bigmatrix()`: Stream variable
  correlations into a destination big.matrix.

- `pca_variable_contributions_stream_bigmatrix()`: Stream variable
  contributions into a destination big.matrix.

## Examples

``` r
set.seed(456)
mat <- bigmemory::as.big.matrix(matrix(rnorm(30), nrow = 6))
ncomp <- 2
rotation_store <- bigmemory::big.matrix(ncol(mat), ncomp, type = "double")
pca_stream <- pca_stream_bigmatrix(mat, xpRotation = rotation_store, ncomp = ncomp)
score_store <- bigmemory::big.matrix(nrow(mat), ncomp, type = "double")
pca_scores_stream_bigmatrix(
    mat,
    score_store,
    pca_stream$rotation,
    pca_stream$center,
    pca_stream$scale,
    ncomp = ncomp
)
#> <pointer: 0x1290c8a40>
loadings_store <- bigmemory::big.matrix(ncol(mat), ncomp, type = "double")
pca_variable_loadings_stream_bigmatrix(
    pca_stream$rotation_stream_bigmatrix,
    pca_stream$sdev,
    loadings_store
)
#> <pointer: 0x109ff98c0>
correlation_store <- bigmemory::big.matrix(ncol(mat), ncomp, type = "double")
pca_variable_correlations_stream_bigmatrix(
    pca_stream$rotation_stream_bigmatrix,
    pca_stream$sdev,
    pca_stream$column_sd,
    pca_stream$scale,
    correlation_store
)
#> <pointer: 0x109fa1190>
contribution_store <- bigmemory::big.matrix(ncol(mat), ncomp, type = "double")
pca_variable_contributions_stream_bigmatrix(
    loadings_store,
    contribution_store
)
#> <pointer: 0x109f8ff40>
```
