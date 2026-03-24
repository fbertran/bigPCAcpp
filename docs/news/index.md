# Changelog

## bigPCAcpp 0.9.1

- Added irlba to benchmarks.

## bigPCAcpp 0.9.0

CRAN release: 2025-10-20

- Added vignettes, readme.

- Initial CRAN submission.

## bigPCAcpp 0.5.0

- Renamed the R implementation of scalable PCA to
  [`pca_spca_R()`](https://fbertran.github.io/bigPCAcpp/reference/pca_spca.md)
  and added a C++-backed
  [`pca_spca()`](https://fbertran.github.io/bigPCAcpp/reference/pca_spca.md)
  along with
  [`pca_spca_stream_bigmatrix()`](https://fbertran.github.io/bigPCAcpp/reference/pca_stream_bigmatrix.md)
  for
  [`bigmemory::big.matrix`](https://rdrr.io/pkg/bigmemory/man/big.matrix.html)
  workflows, retaining the block power iteration algorithm of Elgamal et
  al. (2015).

- Added
  [`pca_spca()`](https://fbertran.github.io/bigPCAcpp/reference/pca_spca.md)
  implementing the scalable PCA algorithm of Elgamal et al. (2015) with
  streaming block power iterations for large matrices.

## bigPCAcpp 0.1.0

- Added an iteratively reweighted singular value decomposition backend
  for
  [`pca_robust()`](https://fbertran.github.io/bigPCAcpp/reference/pca_robust.md)
  that exposes the final row weights and iteration count.

- Implemented a compiled
  [`svd_robust()`](https://fbertran.github.io/bigPCAcpp/reference/svd_robust.md)
  helper with an
  [`svd_robust_R()`](https://fbertran.github.io/bigPCAcpp/reference/svd_robust_R.md)
  reference implementation.
