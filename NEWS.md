# bigPCAcpp 0.9.2

* Added an optional `filematrix` row-block provider for benchmark-side
  validation artifacts without changing the existing `bigmemory` PCA backends.

* Added true `filematrix` support for streaming SPCA via
  `pca_spca_stream_filematrix()` and dense block-wise score projection via
  `pca_scores_stream_filematrix()`.

* Added exact `filematrix` PCA via `pca_stream_filematrix()` for moderate-p
  benchmark validation. This backend streams row blocks but forms a `p x p`
  covariance matrix, so it is not intended for very-wide p; use
  `pca_spca_stream_filematrix()` for very-wide `filematrix` workflows.

# bigPCAcpp 0.9.1

* Added irlba to benchmarks.

# bigPCAcpp 0.9.0

* Added vignettes, readme.

* Initial CRAN submission.

# bigPCAcpp 0.5.0

* Renamed the R implementation of scalable PCA to `pca_spca_R()` and added a
  C++-backed `pca_spca()` along with `pca_spca_stream_bigmatrix()` for
  `bigmemory::big.matrix` workflows, retaining the block power iteration
  algorithm of Elgamal et al. (2015).

* Added `pca_spca()` implementing the scalable PCA algorithm of Elgamal et 
  al. (2015) with streaming block power iterations for large matrices.

# bigPCAcpp 0.1.0

* Added an iteratively reweighted singular value decomposition backend for
  `pca_robust()` that exposes the final row weights and iteration count.

* Implemented a compiled `svd_robust()` helper with an `svd_robust_R()`
  reference implementation.
