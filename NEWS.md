# bigPCA 0.1.0

* Initial CRAN submission.

* Renamed the R implementation of scalable PCA to `pca_spca_R()` and added a
  C++-backed `pca_spca()` along with `pca_spca_stream_bigmatrix()` for
  `bigmemory::big.matrix` workflows, retaining the block power iteration
  algorithm of Elgamal et al. (2015).

* Added `pca_spca()` implementing the scalable PCA algorithm of Elgamal et al.
  (2015) with streaming block power iterations for large matrices.

* Added an iteratively reweighted singular value decomposition backend for
  `pca_robust()` that exposes the final row weights and iteration count.

* Implemented a compiled `svd_robust()` helper with an `svd_robust_R()`
  reference implementation.
