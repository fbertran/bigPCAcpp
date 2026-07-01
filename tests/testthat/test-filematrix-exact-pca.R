align_columns_filematrix_exact <- function(reference, target) {
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

make_exact_test_filematrix <- function(x, label) {
  base <- file.path(
    tempdir(),
    paste("bigpcacpp_exact", label, Sys.getpid(), sample.int(100000000L, 1L), sep = "_")
  )
  fm <- filematrix::fm.create(filenamebase = base, nrow = nrow(x), ncol = ncol(x), type = "double")
  fm[,] <- x
  fm
}

test_that("pca_stream_filematrix matches prcomp with centering", {
  skip_if_not_installed("filematrix")

  set.seed(1101)
  dense <- matrix(rnorm(80), nrow = 20, ncol = 4)
  colnames(dense) <- paste0("v", seq_len(ncol(dense)))
  rownames(dense) <- paste0("r", seq_len(nrow(dense)))

  fm <- make_exact_test_filematrix(dense, "centered")
  on.exit(filematrix::closeAndDeleteFiles(fm), add = TRUE)
  colnames(fm) <- colnames(dense)
  rownames(fm) <- rownames(dense)

  res <- pca_stream_filematrix(
    fm,
    center = TRUE,
    scale = FALSE,
    ncomp = 3L,
    chunk_size = 6L,
    return_scores = TRUE
  )
  pr <- prcomp(dense, center = TRUE, scale. = FALSE)

  expect_s3_class(res, c("bigpca_stream_filematrix", "bigpca"))
  expect_equal(attr(res, "backend"), "stream_filematrix")
  expect_equal(res$nobs, nrow(dense))
  expect_equal(res$center, pr$center, tolerance = 1e-10)
  expect_null(res$scale)
  expect_equal(res$column_sd, apply(dense, 2, stats::sd), tolerance = 1e-10)
  expect_equal(unname(res$sdev), unname(pr$sdev[seq_len(3)]), tolerance = 1e-8)
  expect_equal(res$covariance, stats::cov(dense), tolerance = 1e-10)

  aligned_rotation <- align_columns_filematrix_exact(pr$rotation[, seq_len(3)], res$rotation)
  expect_equal(aligned_rotation, pr$rotation[, seq_len(3)], tolerance = 1e-8)

  aligned_scores <- align_columns_filematrix_exact(pr$x[, seq_len(3)], res$scores)
  expect_equal(aligned_scores, pr$x[, seq_len(3)], tolerance = 1e-8)
})

test_that("pca_stream_filematrix matches prcomp with centering and scaling", {
  skip_if_not_installed("filematrix")

  set.seed(1102)
  dense <- matrix(rnorm(150), nrow = 25, ncol = 6)
  fm <- make_exact_test_filematrix(dense, "scaled")
  on.exit(filematrix::closeAndDeleteFiles(fm), add = TRUE)

  res <- pca_stream_filematrix(
    fm,
    center = TRUE,
    scale = TRUE,
    ncomp = 4L,
    chunk_size = 7L,
    return_scores = TRUE
  )
  pr <- prcomp(dense, center = TRUE, scale. = TRUE)

  expect_equal(res$center, pr$center, tolerance = 1e-10)
  expect_equal(res$scale, pr$scale, tolerance = 1e-10)
  expect_equal(unname(res$sdev), unname(pr$sdev[seq_len(4)]), tolerance = 1e-8)

  scaled_dense <- scale(dense, center = pr$center, scale = pr$scale)
  expect_equal(res$covariance, stats::cov(scaled_dense), tolerance = 1e-10)

  aligned_rotation <- align_columns_filematrix_exact(pr$rotation[, seq_len(4)], res$rotation)
  expect_equal(aligned_rotation, pr$rotation[, seq_len(4)], tolerance = 1e-8)
  aligned_scores <- align_columns_filematrix_exact(pr$x[, seq_len(4)], res$scores)
  expect_equal(aligned_scores, pr$x[, seq_len(4)], tolerance = 1e-8)
})

test_that("pca_scores_stream_filematrix reuses exact filematrix rotation", {
  skip_if_not_installed("filematrix")

  set.seed(1103)
  dense <- matrix(rnorm(90), nrow = 18, ncol = 5)
  fm <- make_exact_test_filematrix(dense, "scores")
  on.exit(filematrix::closeAndDeleteFiles(fm), add = TRUE)

  res <- pca_stream_filematrix(
    fm,
    center = TRUE,
    scale = TRUE,
    ncomp = 3L,
    chunk_size = 5L,
    return_scores = TRUE
  )
  scores <- pca_scores_stream_filematrix(
    fm,
    rotation = res$rotation,
    center = res$center,
    scale = res$scale,
    ncomp = 2L,
    chunk_size = 4L
  )
  expected <- scale(dense, center = res$center, scale = res$scale) %*% res$rotation[, seq_len(2), drop = FALSE]

  expect_equal(scores, expected, tolerance = 1e-10)
  expect_equal(res$scores[, seq_len(2), drop = FALSE], expected, tolerance = 1e-10)
})

test_that("pca_stream_filematrix handles a moderate-p benchmark case", {
  skip_if_not_installed("filematrix")

  set.seed(1104)
  dense <- matrix(rnorm(40 * 30), nrow = 40, ncol = 30)
  fm <- make_exact_test_filematrix(dense, "moderate_p")
  on.exit(filematrix::closeAndDeleteFiles(fm), add = TRUE)

  res <- pca_stream_filematrix(
    fm,
    center = TRUE,
    scale = FALSE,
    ncomp = 5L,
    chunk_size = 9L
  )
  pr <- prcomp(dense, center = TRUE, scale. = FALSE)

  expect_s3_class(res, c("bigpca_stream_filematrix", "bigpca"))
  expect_equal(unname(res$sdev), unname(pr$sdev[seq_len(5)]), tolerance = 1e-8)
  aligned_rotation <- align_columns_filematrix_exact(pr$rotation[, seq_len(5)], res$rotation)
  expect_equal(aligned_rotation, pr$rotation[, seq_len(5)], tolerance = 1e-8)
})

test_that("pca_stream_filematrix validates exact covariance safeguards", {
  skip_if_not_installed("filematrix")

  dense <- matrix(rnorm(24), nrow = 6, ncol = 4)
  fm <- make_exact_test_filematrix(dense, "guard")
  on.exit(filematrix::closeAndDeleteFiles(fm), add = TRUE)

  expect_error(
    pca_stream_filematrix(fm, center = FALSE, scale = TRUE),
    "Scaling requires centring"
  )
  expect_error(
    pca_stream_filematrix(fm, ncomp = NA_integer_),
    "`ncomp` must select at least one component",
    fixed = TRUE
  )

  withr::local_options(list(bigPCAcpp.filematrix.max_cov_gb = 0))
  expect_error(
    pca_stream_filematrix(fm),
    "Refusing exact filematrix PCA"
  )
})

test_that("pca_stream_filematrix rejects non-finite filematrix blocks", {
  skip_if_not_installed("filematrix")

  dense <- matrix(rnorm(24), nrow = 6, ncol = 4)
  dense[2, 3] <- Inf
  fm <- make_exact_test_filematrix(dense, "nonfinite")
  on.exit(filematrix::closeAndDeleteFiles(fm), add = TRUE)

  expect_error(
    pca_stream_filematrix(fm, chunk_size = 2L),
    "Provider block contains non-finite values"
  )
})
