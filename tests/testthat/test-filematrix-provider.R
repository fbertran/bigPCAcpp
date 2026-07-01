test_that("filematrix row provider reads deterministic row and column blocks", {
  skip_if_not_installed("filematrix")

  base <- file.path(
    tempdir(),
    paste("bigpcacpp_filematrix_provider", Sys.getpid(), sample.int(100000000L, 1L), sep = "_")
  )

  X <- matrix(seq_len(42), nrow = 6, ncol = 7)
  storage.mode(X) <- "double"

  fm <- filematrix::fm.create(filenamebase = base, nrow = nrow(X), ncol = ncol(X), type = "double")
  on.exit(filematrix::closeAndDeleteFiles(fm), add = TRUE)
  fm[,] <- X

  expect_true(is_filematrix_object(fm))
  expect_false(is_filematrix_object(X))

  provider <- make_filematrix_row_provider(fm, chunk_size = 2L)

  expect_s3_class(provider, "filematrix_row_provider")
  expect_equal(provider$storage_type, "filematrix")
  expect_equal(provider$chunk_size, 2L)
  expect_equal(provider$nrow(), nrow(X))
  expect_equal(provider$ncol(), ncol(X))

  expect_equal(provider$get_rows(1L, 2L), X[1:2, , drop = FALSE])
  expect_equal(provider$get_rows(3L, 4L), X[3:4, , drop = FALSE])
  expect_equal(provider$get_rows(6L, 6L), X[6L, , drop = FALSE])

  expect_equal(provider$get_cols(2L, 4L), X[, 2:4, drop = FALSE])
  expect_equal(provider$get_block(2L, 5L, 3L, 6L), X[2:5, 3:6, drop = FALSE])
})

test_that("filematrix provider returns double matrices for integer storage", {
  skip_if_not_installed("filematrix")

  base <- file.path(
    tempdir(),
    paste("bigpcacpp_filematrix_provider_integer", Sys.getpid(), sample.int(100000000L, 1L), sep = "_")
  )

  X <- matrix(seq_len(12), nrow = 4, ncol = 3)
  fm <- filematrix::fm.create(filenamebase = base, nrow = nrow(X), ncol = ncol(X), type = "integer")
  on.exit(filematrix::closeAndDeleteFiles(fm), add = TRUE)
  fm[,] <- X

  provider <- make_filematrix_row_provider(fm, chunk_size = 2L)
  block <- provider$get_rows(2L, 3L)

  expect_equal(block, X[2:3, , drop = FALSE])
  expect_true(is.double(block))
})

test_that("filematrix provider rejects unsupported non-numeric storage", {
  skip_if_not_installed("filematrix")

  base <- file.path(
    tempdir(),
    paste("bigpcacpp_filematrix_provider_logical", Sys.getpid(), sample.int(100000000L, 1L), sep = "_")
  )

  fm <- filematrix::fm.create(filenamebase = base, nrow = 2, ncol = 2, type = "logical")
  on.exit(filematrix::closeAndDeleteFiles(fm), add = TRUE)
  fm[,] <- matrix(c(TRUE, FALSE, TRUE, FALSE), nrow = 2)

  expect_error(make_filematrix_row_provider(fm), "Only numeric filematrix storage is supported")
})

test_that("filematrix provider validates inputs and bounds", {
  skip_if_not_installed("filematrix")

  expect_error(make_filematrix_row_provider(matrix(1, 2, 2)), "filematrix object")

  base <- file.path(
    tempdir(),
    paste("bigpcacpp_filematrix_provider_bounds", Sys.getpid(), sample.int(100000000L, 1L), sep = "_")
  )

  fm <- filematrix::fm.create(filenamebase = base, nrow = 3, ncol = 2, type = "double")
  on.exit(filematrix::closeAndDeleteFiles(fm), add = TRUE)
  fm[,] <- matrix(rnorm(6), 3, 2)

  provider <- make_filematrix_row_provider(fm)

  expect_error(provider$get_rows(0L, 1L), "positive scalar integer")
  expect_error(provider$get_rows(2L, 1L), "start must be <=")
  expect_error(provider$get_rows(1L, 4L), "outside provider dimensions")
  expect_error(provider$get_block(1L, 2L, 1L, 3L), "outside provider dimensions")
})

test_that("filematrix provider is benchmark groundwork, not a PCA backend", {
  skip_if_not_installed("filematrix")

  expect_true("make_filematrix_row_provider" %in% getNamespaceExports("bigPCAcpp"))
  expect_true("pca_spca_stream_filematrix" %in% getNamespaceExports("bigPCAcpp"))
  expect_true("pca_scores_stream_filematrix" %in% getNamespaceExports("bigPCAcpp"))
})
