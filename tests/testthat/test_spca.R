skip_if_not_installed("stats")

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
  if (!is.null(rownames(ref_mat)) && is.null(rownames(tgt_mat))) {
    rownames(tgt_mat) <- rownames(ref_mat)
  }
  if (!is.null(colnames(ref_mat))) {
    tgt_names <- colnames(tgt_mat)
    if (is.null(tgt_names)) {
      tgt_names <- rep_len("", ncol(tgt_mat))
    }
    tgt_names[seq_len(cols)] <- colnames(ref_mat)[seq_len(cols)]
    colnames(tgt_mat) <- tgt_names
  }
  tgt_mat
}

test_that("pca_spca approximates dense PCA", {
  set.seed(2025)
  dense <- matrix(rnorm(400), nrow = 100, ncol = 4)
  pr <- prcomp(dense, center = TRUE, scale. = FALSE)

  spca <- pca_spca_R(
    dense,
    ncomp = 3,
    center = TRUE,
    scale = FALSE,
    block_size = 32L,
    max_iter = 80L,
    tol = 1e-4,
    seed = 123,
    return_scores = TRUE
  )

  expect_s3_class(spca, c("bigpca_spca_r", "bigpca"))
  expect_equal(attr(spca, "backend"), "spca_r")
  expect_equal(spca$nobs, nrow(dense))
  expect_null(spca$scale)
  expect_equal(spca$center, colMeans(dense), tolerance = 1e-8)
  expect_equal(spca$column_sd, apply(dense, 2, stats::sd), tolerance = 1e-6)
  expect_equal(unname(spca$sdev), unname(pr$sdev[seq_len(3)]), tolerance = 1e-5)

  rotation_ref <- pr$rotation[, seq_len(3), drop = FALSE]
  rotation_mat <- align_columns(rotation_ref, spca$rotation)
  expect_equal(rotation_mat, rotation_ref, tolerance = 1e-3)

  scores_ref <- pr$x[, seq_len(3), drop = FALSE]
  aligned_scores <- align_columns(scores_ref, spca$scores)
  expect_equal(aligned_scores, scores_ref, tolerance = 1e-3)

  total_var <- sum(apply(dense, 2, stats::sd)^2)
  expect_equal(spca$explained_variance, spca$eigenvalues / total_var, tolerance = 1e-5)
  expect_true(isTRUE(attr(spca, "converged")))
  expect_gt(attr(spca, "iterations"), 0)
})

test_that("pca_spca supports scaling to unit variance", {
  set.seed(42)
  dense <- matrix(rnorm(150), nrow = 50, ncol = 3)
  pr <- prcomp(dense, center = TRUE, scale. = TRUE)

  spca <- pca_spca_R(
    dense,
    ncomp = 2,
    center = TRUE,
    scale = TRUE,
    block_size = 16L,
    max_iter = 80L,
    tol = 1e-4,
    seed = 99,
    return_scores = TRUE
  )

  centered <- sweep(dense, 2, colMeans(dense), "-")
  sd_vec <- apply(centered, 2, stats::sd)
  expect_equal(spca$scale, sd_vec, tolerance = 1e-6)
  expect_equal(spca$column_sd, sd_vec, tolerance = 1e-6)
  expect_equal(unname(spca$sdev), unname(pr$sdev[seq_len(2)]), tolerance = 1e-5)

  rotation_ref <- pr$rotation[, seq_len(2), drop = FALSE]
  rotation_mat <- align_columns(rotation_ref, spca$rotation)
  expect_equal(rotation_mat, rotation_ref, tolerance = 2e-3)

  scores_ref <- pr$x[, seq_len(2), drop = FALSE]
  aligned_scores <- align_columns(scores_ref, spca$scores)
  expect_equal(aligned_scores, scores_ref, tolerance = 2e-3)
  expect_true(isTRUE(attr(spca, "converged")))
})

test_that("pca_spca processes big.matrix inputs via C++ backend", {
  skip_if_not_installed("bigmemory")
  set.seed(303)
  dense <- matrix(rnorm(200), nrow = 50, ncol = 4)
  bm <- bigmemory::as.big.matrix(dense)

  cpp_res <- pca_spca(
    bm,
    ncomp = 3,
    center = TRUE,
    scale = FALSE,
    block_size = 32L,
    max_iter = 60L,
    tol = 1e-4,
    seed = 7L,
    return_scores = TRUE
  )

  ref_res <- pca_spca_R(
    dense,
    ncomp = 3,
    center = TRUE,
    scale = FALSE,
    block_size = 32L,
    max_iter = 60L,
    tol = 1e-4,
    seed = 7L,
    return_scores = TRUE
  )

  expect_s3_class(cpp_res, c("bigpca_spca_bigmemory", "bigpca"))
  expect_equal(attr(cpp_res, "backend"), "spca_bigmemory")
  expect_equal(cpp_res$nobs, nrow(dense))
  expect_equal(cpp_res$center, ref_res$center, tolerance = 1e-8)
  expect_equal(unname(cpp_res$sdev), unname(ref_res$sdev), tolerance = 1e-5)

  aligned_rotation <- align_columns(ref_res$rotation, cpp_res$rotation)
  expect_equal(aligned_rotation, ref_res$rotation, tolerance = 1e-3)

  aligned_scores <- align_columns(ref_res$scores, cpp_res$scores)
  expect_equal(aligned_scores, ref_res$scores, tolerance = 1e-3)
  expect_true(isTRUE(attr(cpp_res, "converged")))
})

test_that("pca_spca_stream_bigmatrix writes rotation into big.matrix", {
  skip_if_not_installed("bigmemory")
  set.seed(404)
  dense <- matrix(rnorm(300), nrow = 60, ncol = 5)
  bm <- bigmemory::as.big.matrix(dense)
  rot_store <- bigmemory::big.matrix(nrow = ncol(dense), ncol = 2, type = "double")

  streamed <- pca_spca_stream_bigmatrix(
    bm,
    xpRotation = rot_store@address,
    ncomp = 2,
    center = TRUE,
    scale = FALSE,
    block_size = 32L,
    max_iter = 40L,
    tol = 1e-4,
    seed = 11L,
    return_scores = FALSE
  )

  expect_s3_class(streamed, c("bigpca_spca_stream_bigmatrix", "bigpca"))
  expect_equal(attr(streamed, "backend"), "spca_stream_bigmatrix")
  expect_true("rotation_stream_bigmatrix" %in% names(streamed))
  stored_rotation <- bigmemory::as.matrix(rot_store)
  aligned_rotation <- align_columns(streamed$rotation, stored_rotation)
  expect_equal(aligned_rotation, streamed$rotation, tolerance = 1e-8)
})
