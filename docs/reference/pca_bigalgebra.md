# bigalgebra-backed PCA helpers

Variants of the PCA helpers that stream results directly into
[`bigmemory::big.matrix`](https://rdrr.io/pkg/bigmemory/man/big.matrix.html)
objects, enabling file-backed workflows without materialising dense R
matrices.

## Usage

``` r
pca_bigalgebra(
  xpMat,
  xpRotation = NULL,
  center = TRUE,
  scale = FALSE,
  ncomp = -1L,
  block_size = 1024L
)

pca_scores_bigalgebra(
  xpMat,
  xpDest,
  rotation,
  center,
  scale,
  ncomp = -1L,
  block_size = 1024L
)

pca_variable_loadings_bigalgebra(xpRotation, sdev, xpDest)
```

## Arguments

- xpRotation:

  Optionally, either a
  [`bigmemory::big.matrix`](https://rdrr.io/pkg/bigmemory/man/big.matrix.html)
  or external pointer referencing a destination big.matrix that receives
  the rotation matrix.

- xpDest:

  Either a `big.matrix` or external pointer referencing the destination
  `big.matrix` that stores the computed quantity.

- xpLoadings:

  For `pca_variable_contributions_bigalgebra()`, the loadings matrix
  supplied as a `big.matrix` or external pointer.

## Value

For `pca_bigalgebra()`, the same `bigpca` object as `pca_bigmatrix()`
with the addition of a `rotation_bigalgebra` element referencing the
populated `big.matrix` when `xpRotation` is supplied.

The external pointer supplied in `xpDest`, invisibly.

## Functions

- `pca_scores_bigalgebra()`: Stream PCA scores into a destination
  big.matrix.

- `pca_variable_loadings_bigalgebra()`: Populate big.matrix objects with
  derived variable diagnostics.

## Examples

``` r
set.seed(456)
mat <- bigmemory::as.big.matrix(matrix(rnorm(30), nrow = 6))
rotation_store <- bigmemory::big.matrix(ncol(mat), ncol(mat), type = "double")
pca_stream <- pca_bigalgebra(mat, xpRotation = rotation_store, ncomp = 2)
#> Error: rotation big.matrix has incompatible dimensions
score_store <- bigmemory::big.matrix(nrow(mat), 2, type = "double")
pca_scores_bigalgebra(mat, score_store, pca_stream$rotation, pca_stream$center,
    pca_stream$scale, ncomp = 2)
#> Error: object 'pca_stream' not found
```
