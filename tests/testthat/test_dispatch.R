skip_if_not_installed("bigmemory")

create_double_big_matrix <- function(mat) {
  bm <- bigmemory::big.matrix(
    nrow = nrow(mat),
    ncol = ncol(mat),
    type = "double"
  )
  bm[,] <- mat
  bm
}

maybe_vec <- function(x) {
  if (is.null(x)) numeric() else x
}

align_columns <- function(reference, target) {
  ref_mat <- as.matrix(reference)
  tgt_mat <- as.matrix(target)
  cols <- min(ncol(ref_mat), ncol(tgt_mat))
  if (cols == 0) {
    return(tgt_mat)
  }
  for (j in seq_len(cols)) {
    ref <- ref_mat[, j]
    tgt <- tgt_mat[, j]
    if (sum((ref - tgt)^2) > sum((ref + tgt)^2)) {
      tgt_mat[, j] <- -tgt
    }
  }
  if (!is.null(colnames(ref_mat))) {
    colnames(tgt_mat) <- colnames(ref_mat)
  }
  if (!is.null(rownames(ref_mat))) {
    rownames(tgt_mat) <- rownames(ref_mat)
  }
  tgt_mat
}

test_that("double big.matrix wrappers dispatch to double implementations", {
  set.seed(2025)
  dense <- matrix(rnorm(36, mean = 1e8, sd = 1e4), nrow = 12, ncol = 3)

  bm_double <- create_double_big_matrix(dense)

  base_pca <- prcomp(dense, center = TRUE, scale. = FALSE)
  pca_res <- pca_bigmatrix(bm_double, center = TRUE, scale = FALSE, ncomp = 3L)
  expect_equal(pca_res$sdev, base_pca$sdev, tolerance = 1e-8)
  aligned_rotation <- align_columns(base_pca$rotation, pca_res$rotation[])
  expect_equal(aligned_rotation, base_pca$rotation, tolerance = 1e-8)

  spca_res <- pca_spca(
    bm_double,
    ncomp = 3L,
    center = TRUE,
    scale = FALSE,
    block_size = 4L,
    max_iter = 200L,
    tol = 1e-7,
    return_scores = TRUE,
    verbose = FALSE
  )
  expect_equal(unname(spca_res$sdev), unname(pca_res$sdev), tolerance = 1e-6)
  aligned_spca <- align_columns(pca_res$rotation, spca_res$rotation[])
  expect_equal(unname(aligned_spca), unname(pca_res$rotation[]), tolerance = 1e-6)

  scores_res <- pca_scores_bigmatrix(
    bm_double,
    pca_res$rotation,
    maybe_vec(pca_res$center),
    maybe_vec(pca_res$scale),
    ncomp = 3L,
    block_size = 4L
  )
  aligned_scores <- align_columns(base_pca$x, scores_res[])
  expect_equal(aligned_scores, base_pca$x, tolerance = 1e-8)

  svd_res <- svd_bigmatrix(bm_double, nu = 3L, nv = 3L, block_size = 4L)
  base_svd <- svd(dense)
  expect_equal(svd_res$d, base_svd$d, tolerance = 1e-8)
  aligned_u <- align_columns(base_svd$u, svd_res$u)
  aligned_v <- align_columns(base_svd$v, svd_res$v)
  expect_equal(aligned_u, base_svd$u, tolerance = 1e-8)
  expect_equal(aligned_v, base_svd$v, tolerance = 1e-8)
})
