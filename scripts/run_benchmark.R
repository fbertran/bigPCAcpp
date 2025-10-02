suppressPackageStartupMessages({
  library(bigmemory)
  if (requireNamespace("bigPCAcpp", quietly = TRUE)) {
    library(bigPCAcpp)
  } else {
    if (!requireNamespace("pkgload", quietly = TRUE)) {
      stop("bigPCAcpp must be installed or pkgload must be available", call. = FALSE)
    }
    pkgload::load_all(".")
  }
})

sizes <- list(
  small = list(rows = 1000L, cols = 50L),
  medium = list(rows = 5000L, cols = 100L),
  large = list(rows = 20000L, cols = 200L),
  xlarge = list(rows = 50000L, cols = 300L),
  xxlarge = list(rows = 100000L, cols = 500L),
  xxxlarge = list(rows = 100000L, cols = 2000L)
)

method_runners <- list(
  classical = function(mats, ncomp) {
    pca_bigmatrix(mats$big, center = TRUE, scale = TRUE, ncomp = ncomp)
  },
  streaming = function(mats, ncomp) {
    pca_stream_bigmatrix(mats$big, center = TRUE, scale = TRUE, ncomp = ncomp)
  },
  scalable = function(mats, ncomp) {
    pca_spca(
      mats$big,
      ncomp = ncomp,
      center = TRUE,
      scale = TRUE,
      block_size = 2048L,
      max_iter = 25L,
      tol = 1e-4,
      seed = 42L,
      return_scores = FALSE,
      verbose = FALSE
    )
  },
  prcomp = function(mats, ncomp) {
    stats::prcomp(
      mats$dense,
      center = TRUE,
      scale. = TRUE,
      rank. = ncomp
    )
  }
)

replicates_for <- function(rows) {
  if (rows <= 5000L) {
    20L
  } else if (rows <= 20000L) {
    20L
  } else {
    10L
  }
}

results <- list()
row_id <- 1L
set.seed(123)

for (dataset_name in names(sizes)) {
  dims <- sizes[[dataset_name]]
  message(sprintf("Generating dataset '%s' with %d rows and %d columns", dataset_name, dims$rows, dims$cols))
  mat <- matrix(rnorm(dims$rows * dims$cols), nrow = dims$rows, ncol = dims$cols)
  big_mat <- bigmemory::as.big.matrix(mat, type = "double")
  ncomp <- min(10L, dims$cols)
  reps <- replicates_for(dims$rows)
  inputs <- list(dense = mat, big = big_mat)
  
  for (method_name in names(method_runners)) {
    runner <- method_runners[[method_name]]
    for (rep in seq_len(reps)) {
      gc()
      gc()
      message(sprintf("Running %s (replicate %d/%d) on %s", method_name, rep, reps, dataset_name))
      res <- NULL
      timing <- system.time({
        res <<- tryCatch(
          runner(inputs, ncomp),
          error = function(e) e
        )
      })
      success <- !inherits(res, "error")
      backend <- if (success) {
        backend_val <- attr(res, "backend", exact = TRUE)
        if (is.null(backend_val) && !is.null(res$backend)) {
          res$backend
        } else {
          backend_val
        }
      } else {
        NA_character_
      }
      iterations <- if (success) {
        iter <- attr(res, "iterations", exact = TRUE)
        if (is.null(iter)) NA_integer_ else as.integer(iter)
      } else {
        NA_integer_
      }
      converged <- if (success) {
        conv <- attr(res, "converged", exact = TRUE)
        if (is.null(conv)) NA else as.logical(conv)
      } else {
        NA
      }
      results[[row_id]] <- data.frame(
        dataset = dataset_name,
        rows = dims$rows,
        cols = dims$cols,
        ncomp = ncomp,
        method = method_name,
        replicate = rep,
        user_time = unname(timing[["user.self"]]),
        system_time = unname(timing[["sys.self"]]),
        elapsed = unname(timing[["elapsed"]]),
        success = success,
        backend = if (is.null(backend)) NA_character_ else as.character(backend),
        iterations = iterations,
        converged = converged,
        error = if (success) NA_character_ else conditionMessage(res),
        stringsAsFactors = FALSE
      )
      row_id <- row_id + 1L
    }
  }
  rm(mat, big_mat)
  gc()
  gc()
}

benchmark_results <- do.call(rbind, results)

if (!dir.exists("data")) {
  dir.create("data")
}

save(benchmark_results, file = file.path("data", "benchmark_results.rda"), compress = "bzip2")
