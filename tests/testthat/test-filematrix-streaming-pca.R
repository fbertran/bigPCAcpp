align_columns_filematrix <- function(reference, target) {
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
  tgt_mat
}

make_test_filematrix <- function(x, label) {
  base <- file.path(
    tempdir(),
    paste("bigpcacpp", label, Sys.getpid(), sample.int(100000000L, 1L), sep = "_")
  )
  fm <- filematrix::fm.create(filenamebase = base, nrow = nrow(x), ncol = ncol(x), type = "double")
  fm[,] <- x
  fm
}

test_that("pca_spca_stream_filematrix matches dense SPCA on small data", {
  skip_if_not_installed("filematrix")

  set.seed(707)
  dense <- matrix(rnorm(240), nrow = 40, ncol = 6)
  colnames(dense) <- paste0("v", seq_len(ncol(dense)))
  rownames(dense) <- paste0("r", seq_len(nrow(dense)))

  fm <- make_test_filematrix(dense, "spca_stream")
  on.exit(filematrix::closeAndDeleteFiles(fm), add = TRUE)
  colnames(fm) <- colnames(dense)
  rownames(fm) <- rownames(dense)

  file_res <- pca_spca_stream_filematrix(
    fm,
    ncomp = 2,
    center = TRUE,
    scale = FALSE,
    chunk_size = 9L,
    max_iter = 80L,
    tol = 1e-5,
    seed = 11L,
    return_scores = TRUE
  )
  dense_res <- pca_spca_R(
    dense,
    ncomp = 2,
    center = TRUE,
    scale = FALSE,
    block_size = 9L,
    max_iter = 80L,
    tol = 1e-5,
    seed = 11L,
    return_scores = TRUE
  )

  expect_s3_class(file_res, c("bigpca_spca_stream_filematrix", "bigpca"))
  expect_equal(attr(file_res, "backend"), "spca_stream_filematrix")
  expect_true(isTRUE(attr(file_res, "experimental")))
  expect_equal(attr(file_res, "storage_type"), "filematrix")
  expect_null(file_res$covariance)
  expect_true(all(is.finite(file_res$rotation)))
  expect_true(all(is.finite(file_res$sdev)))
  expect_equal(file_res$center, dense_res$center, tolerance = 1e-10)
  expect_equal(file_res$column_sd, dense_res$column_sd, tolerance = 1e-10)
  expect_equal(unname(file_res$sdev), unname(dense_res$sdev), tolerance = 1e-6)

  aligned_rotation <- align_columns_filematrix(dense_res$rotation, file_res$rotation)
  expect_equal(aligned_rotation, dense_res$rotation, tolerance = 1e-6)
  aligned_scores <- align_columns_filematrix(dense_res$scores, file_res$scores)
  expect_equal(aligned_scores, dense_res$scores, tolerance = 1e-6)
})

test_that("pca_spca dispatches filematrix inputs to streaming backend", {
  skip_if_not_installed("filematrix")

  dense <- matrix(rnorm(120), nrow = 30, ncol = 4)
  fm <- make_test_filematrix(dense, "spca_dispatch")
  on.exit(filematrix::closeAndDeleteFiles(fm), add = TRUE)

  res <- pca_spca(
    fm,
    ncomp = 2,
    center = TRUE,
    scale = FALSE,
    block_size = 8L,
    max_iter = 30L,
    tol = 1e-4,
    seed = 5L
  )

  expect_s3_class(res, c("bigpca_spca_stream_filematrix", "bigpca"))
  expect_equal(attr(res, "backend"), "spca_stream_filematrix")
})

test_that("pca_scores_stream_filematrix streams scores by row block", {
  skip_if_not_installed("filematrix")

  set.seed(808)
  dense <- matrix(rnorm(150), nrow = 25, ncol = 6)
  fm <- make_test_filematrix(dense, "score_stream")
  on.exit(filematrix::closeAndDeleteFiles(fm), add = TRUE)

  rotation <- qr.Q(qr(matrix(rnorm(ncol(dense) * 3), nrow = ncol(dense), ncol = 3)))
  colnames(rotation) <- paste0("PC", seq_len(ncol(rotation)))
  center <- colMeans(dense)
  scale <- apply(sweep(dense, 2, center, "-"), 2, stats::sd)
  scale[scale == 0] <- 1

  scores <- pca_scores_stream_filematrix(
    fm,
    rotation = rotation,
    center = center,
    scale = scale,
    chunk_size = 7L,
    sink = "r"
  )
  expected <- sweep(sweep(dense, 2, center, "-"), 2, scale, "/") %*% rotation

  expect_equal(scores, expected, tolerance = 1e-10)

  scan_info <- pca_scores_stream_filematrix(
    fm,
    rotation = rotation,
    center = center,
    scale = scale,
    chunk_size = 6L,
    sink = "none"
  )
  expect_equal(scan_info$nobs, nrow(dense))
  expect_equal(scan_info$ncomp, ncol(rotation))
  expect_equal(scan_info$sink, "none")
  expect_equal(scan_info$storage_type, "filematrix")
})

test_that("pca_scores_stream_filematrix can write scores to filematrix", {
  skip_if_not_installed("filematrix")

  dense <- matrix(rnorm(80), nrow = 20, ncol = 4)
  fm <- make_test_filematrix(dense, "score_stream_source")
  on.exit(filematrix::closeAndDeleteFiles(fm), add = TRUE)

  rotation <- qr.Q(qr(matrix(rnorm(ncol(dense) * 2), nrow = ncol(dense), ncol = 2)))
  out_base <- file.path(
    tempdir(),
    paste("bigpcacpp_score_sink", Sys.getpid(), sample.int(100000000L, 1L), sep = "_")
  )
  scores_fm <- pca_scores_stream_filematrix(
    fm,
    rotation = rotation,
    chunk_size = 5L,
    sink = "filematrix",
    filematrix_base = out_base
  )
  on.exit(filematrix::closeAndDeleteFiles(scores_fm), add = TRUE)

  expect_true(is_filematrix_object(scores_fm))
  expect_equal(scores_fm[, ], dense %*% rotation, tolerance = 1e-10)
})
