# Scalable principal component analysis via streaming power iterations

Implements the scalable PCA (sPCA) procedure of Elgamal et al. (2015),
which uses block power iterations to approximate the leading principal
components while streaming the data in manageable chunks. The algorithm
only requires matrix-vector products, allowing large matrices to be
processed without materialising the full cross-product in memory.

Implements the scalable PCA (sPCA) procedure of Elgamal et al. (2015),
which uses block power iterations to approximate the leading principal
components while streaming the data in manageable chunks. The algorithm
only requires matrix-vector products, allowing large matrices to be
processed without materialising the full cross-product in memory.

Implements the scalable PCA (sPCA) procedure of Elgamal et al. (2015),
which uses block power iterations to approximate the leading principal
components while streaming the data in manageable chunks. The algorithm
only requires matrix-vector products, allowing large matrices to be
processed without materialising the full cross-product in memory.

## Usage

``` r
pca_spca(
  x,
  ncomp = NULL,
  center = TRUE,
  scale = FALSE,
  block_size = 2048L,
  max_iter = 50L,
  tol = 1e-04,
  seed = NULL,
  return_scores = FALSE,
  verbose = FALSE
)

pca_spca(
  x,
  ncomp = NULL,
  center = TRUE,
  scale = FALSE,
  block_size = 2048L,
  max_iter = 50L,
  tol = 1e-04,
  seed = NULL,
  return_scores = FALSE,
  verbose = FALSE
)

pca_spca_R(
  x,
  ncomp = NULL,
  center = TRUE,
  scale = FALSE,
  block_size = 2048L,
  max_iter = 50L,
  tol = 1e-04,
  seed = NULL,
  return_scores = FALSE,
  verbose = FALSE
)

pca_spca(
  x,
  ncomp = NULL,
  center = TRUE,
  scale = FALSE,
  block_size = 2048L,
  max_iter = 50L,
  tol = 1e-04,
  seed = NULL,
  return_scores = FALSE,
  verbose = FALSE
)

pca_spca_R(
  x,
  ncomp = NULL,
  center = TRUE,
  scale = FALSE,
  block_size = 2048L,
  max_iter = 50L,
  tol = 1e-04,
  seed = NULL,
  return_scores = FALSE,
  verbose = FALSE
)
```

## Arguments

- x:

  A numeric matrix, data frame,
  [`bigmemory::big.matrix`](https://rdrr.io/pkg/bigmemory/man/big.matrix.html),
  or an external pointer referencing a big.matrix. The input is
  processed in row-wise blocks so that large matrices can be analysed
  without creating dense copies in R memory.

- ncomp:

  Number of principal components to retain. Use `NULL` or a non-positive
  value to keep `min(nrow(x), ncol(x))` components.

- center:

  Logical; should column means be subtracted before performing PCA?

- scale:

  Logical; when `TRUE`, columns are scaled to unit variance after
  centring. Scaling requires `center = TRUE`.

- block_size:

  Number of rows to stream per block when computing column statistics
  and matrix-vector products.

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

## Value

A [`bigpca`](https://fbertran.github.io/bigPCAcpp/reference/bigpca.md)
object containing the approximate PCA solution with the same structure
as
[`pca_bigmatrix()`](https://fbertran.github.io/bigPCAcpp/reference/pca_bigmatrix.md).
The result includes component standard deviations, rotation/loadings,
optional scores, column statistics, and variance summaries. Additional
metadata is stored in `attr(result, "iterations")` (number of iterations
performed), `attr(result, "tolerance")` (requested tolerance), and
`attr(result, "converged")` (logical convergence flag).

A [`bigpca`](https://fbertran.github.io/bigPCAcpp/reference/bigpca.md)
object containing the approximate PCA solution with the same structure
as
[`pca_bigmatrix()`](https://fbertran.github.io/bigPCAcpp/reference/pca_bigmatrix.md).
The result includes component standard deviations, rotation/loadings,
optional scores, column statistics, and variance summaries. Additional
metadata is stored in `attr(result, "iterations")` (number of iterations
performed), `attr(result, "tolerance")` (requested tolerance), and
`attr(result, "converged")` (logical convergence flag).

A [`bigpca`](https://fbertran.github.io/bigPCAcpp/reference/bigpca.md)
object containing the approximate PCA solution with the same structure
as
[`pca_bigmatrix()`](https://fbertran.github.io/bigPCAcpp/reference/pca_bigmatrix.md).
The result includes component standard deviations, rotation/loadings,
optional scores, column statistics, and variance summaries. Additional
metadata is stored in `attr(result, "iterations")` (number of iterations
performed), `attr(result, "tolerance")` (requested tolerance), and
`attr(result, "converged")` (logical convergence flag).

## References

Tarek Elgamal, Maysam Yabandeh, Ashraf Aboulnaga, Waleed Mustafa, and
Mohamed Hefeeda (2015). *sPCA: Scalable Principal Component Analysis for
Big Data on Distributed Platforms*. Proceedings of the 2015 ACM SIGMOD
International Conference on Management of Data.
<doi:10.1145/2723372.2751520>.

Tarek Elgamal, Maysam Yabandeh, Ashraf Aboulnaga, Waleed Mustafa, and
Mohamed Hefeeda (2015). *sPCA: Scalable Principal Component Analysis for
Big Data on Distributed Platforms*. Proceedings of the 2015 ACM SIGMOD
International Conference on Management of Data.
<doi:10.1145/2723372.2751520>.

Tarek Elgamal, Maysam Yabandeh, Ashraf Aboulnaga, Waleed Mustafa, and
Mohamed Hefeeda (2015). *sPCA: Scalable Principal Component Analysis for
Big Data on Distributed Platforms*. Proceedings of the 2015 ACM SIGMOD
International Conference on Management of Data.
<doi:10.1145/2723372.2751520>.
