#' Experimental filematrix PCA helpers
#'
#' @description
#' Experimental streaming helpers for `filematrix` inputs. These functions read
#' the source matrix in row blocks through [make_filematrix_row_provider()],
#' avoiding `bigmemory` memory-mapped access patterns on shared filesystems.
#' `pca_stream_filematrix()` performs exact covariance PCA for moderate-p
#' benchmark validation by forming a `p x p` covariance matrix. Very-wide
#' `filematrix` workflows should use [pca_spca_stream_filematrix()] instead.
#'
#' @inheritParams pca_spca
#' @param X A `filematrix` object.
#' @param center For `pca_stream_filematrix()` and
#'   `pca_spca_stream_filematrix()`, logical; should column means be subtracted
#'   before fitting? For `pca_scores_stream_filematrix()`, `NULL` or a numeric
#'   vector of column centers.
#' @param scale For `pca_stream_filematrix()` and
#'   `pca_spca_stream_filematrix()`, logical; should columns be scaled after
#'   centering? For `pca_scores_stream_filematrix()`, `NULL` or a numeric vector
#'   of column scale factors.
#' @param chunk_size Number of rows to read per block.
#'
#' @return `pca_stream_filematrix()` returns a [`bigpca`] object compatible
#'   with [pca_bigmatrix()] outputs and records the `stream_filematrix` backend.
#'   `pca_spca_stream_filematrix()` returns a [`bigpca`] object compatible with
#'   [pca_spca()] outputs and records the `spca_filematrix` backend.
#'   `pca_scores_stream_filematrix()` returns a dense score matrix computed by
#'   scanning the `filematrix` input in row blocks.
#'
#' @name filematrix_streaming_pca
NULL

#' @rdname filematrix_streaming_pca
#' @export
pca_stream_filematrix <- function(
        X,
        center = TRUE,
        scale = FALSE,
        ncomp = -1L,
        chunk_size = 1024L,
        return_scores = FALSE) {
    if (!requireNamespace("filematrix", quietly = TRUE)) {
        stop("pca_stream_filematrix() requires the optional filematrix package.", call. = FALSE)
    }
    provider <- make_filematrix_row_provider(X, chunk_size = chunk_size)
    .exact_pca_stream_row_provider(
        provider = provider,
        X = X,
        center = center,
        scale = scale,
        ncomp = ncomp,
        return_scores = return_scores
    )
}

#' @rdname filematrix_streaming_pca
#' @export
pca_spca_stream_filematrix <- function(
        X,
        ncomp = NULL,
        center = TRUE,
        scale = FALSE,
        chunk_size = 2048L,
        max_iter = 50L,
        tol = 1e-4,
        seed = NULL,
        return_scores = FALSE,
        verbose = FALSE) {
    if (!requireNamespace("filematrix", quietly = TRUE)) {
        stop("pca_spca_stream_filematrix() requires the optional filematrix package.", call. = FALSE)
    }
    provider <- make_filematrix_row_provider(X, chunk_size = chunk_size)
    run_impl <- function() {
        .spca_stream_row_provider(
            provider = provider,
            X = X,
            ncomp = ncomp,
            center = center,
            scale = scale,
            max_iter = max_iter,
            tol = tol,
            return_scores = return_scores,
            verbose = verbose
        )
    }
    run_with_seed(seed, run_impl)
}

#' @rdname filematrix_streaming_pca
#' @param rotation Numeric rotation matrix whose rows match `ncol(X)`.
#' @param ncomp Number of score columns to return. Use a non-positive value to
#'   keep all columns in `rotation`.
#' @export
pca_scores_stream_filematrix <- function(
        X,
        rotation,
        center = NULL,
        scale = NULL,
        ncomp = -1L,
        chunk_size = 1024L) {
    if (!requireNamespace("filematrix", quietly = TRUE)) {
        stop("pca_scores_stream_filematrix() requires the optional filematrix package.", call. = FALSE)
    }
    provider <- make_filematrix_row_provider(X, chunk_size = chunk_size)
    .scores_stream_row_provider(
        provider = provider,
        X = X,
        rotation = rotation,
        center = center,
        scale = scale,
        ncomp = ncomp
    )
}

.exact_pca_stream_row_provider <- function(provider,
                                           X,
                                           center,
                                           scale,
                                           ncomp,
                                           return_scores) {
    if (!is.logical(center) || length(center) != 1L || is.na(center)) {
        stop("`center` must be TRUE or FALSE", call. = FALSE)
    }
    if (!is.logical(scale) || length(scale) != 1L || is.na(scale)) {
        stop("`scale` must be TRUE or FALSE", call. = FALSE)
    }
    if (scale && !center) {
        stop("Scaling requires centring the variables", call. = FALSE)
    }
    if (!is.logical(return_scores) || length(return_scores) != 1L || is.na(return_scores)) {
        stop("`return_scores` must be TRUE or FALSE", call. = FALSE)
    }

    n <- provider$nrow()
    p <- provider$ncol()
    if (n < 2L) {
        stop("PCA requires at least two observations", call. = FALSE)
    }
    .check_filematrix_covariance_budget(p)

    if (is.null(ncomp)) {
        ncomp <- p
    }
    ncomp <- as.integer(ncomp)
    if (length(ncomp) != 1L) {
        stop("`ncomp` must be a scalar integer", call. = FALSE)
    }
    if (!is.na(ncomp) && ncomp <= 0L) {
        ncomp <- p
    }
    if (is.na(ncomp) || ncomp < 1L) {
        stop("`ncomp` must select at least one component", call. = FALSE)
    }
    if (ncomp > p) {
        ncomp <- p
    }

    block_starts <- seq.int(1L, n, by = provider$chunk_size)
    col_sum <- numeric(p)
    observed <- 0L
    for (start in block_starts) {
        end <- min(start + provider$chunk_size - 1L, n)
        block <- provider$get_rows(start, end)
        .check_finite_numeric_block(block)
        col_sum <- col_sum + colSums(block)
        observed <- observed + nrow(block)
    }
    if (observed != n) {
        stop("Provider returned an unexpected number of rows", call. = FALSE)
    }

    mean_vec <- col_sum / n
    col_ss <- numeric(p)
    for (start in block_starts) {
        end <- min(start + provider$chunk_size - 1L, n)
        block <- provider$get_rows(start, end)
        .check_finite_numeric_block(block)
        centered <- sweep(block, 2, mean_vec, "-", check.margin = FALSE)
        col_ss <- col_ss + colSums(centered^2)
    }

    column_sd <- sqrt(pmax(col_ss / max(1L, n - 1L), 0))
    center_vec <- if (center) mean_vec else NULL
    scale_vec <- if (scale) {
        out <- column_sd
        out[out == 0] <- 1
        out
    } else {
        NULL
    }

    transform_block <- function(block) {
        if (!is.null(center_vec)) {
            block <- sweep(block, 2, center_vec, "-", check.margin = FALSE)
        }
        if (!is.null(scale_vec)) {
            block <- sweep(block, 2, scale_vec, "/", check.margin = FALSE)
        }
        block
    }

    covariance <- matrix(0, nrow = p, ncol = p)
    for (start in block_starts) {
        end <- min(start + provider$chunk_size - 1L, n)
        block <- provider$get_rows(start, end)
        .check_finite_numeric_block(block)
        block <- transform_block(block)
        covariance <- covariance + crossprod(block)
    }
    denom <- n - 1L
    if (denom <= 0L) {
        stop("Unable to compute covariance with fewer than two observations", call. = FALSE)
    }
    covariance <- covariance / denom
    covariance <- (covariance + t(covariance)) / 2

    eig <- eigen(covariance, symmetric = TRUE)
    eigenvalues_all <- pmax(eig$values, 0)
    rotation <- eig$vectors[, seq_len(ncomp), drop = FALSE]
    eigenvalues <- eigenvalues_all[seq_len(ncomp)]
    sdev <- sqrt(eigenvalues)
    total_variance <- sum(eigenvalues_all)
    explained <- if (total_variance > 0) eigenvalues / total_variance else rep(0, length(eigenvalues))
    cumulative <- cumsum(explained)

    colnames(rotation) <- paste0("PC", seq_len(ncomp))
    filematrix_colnames <- .filematrix_colnames(X)
    if (!is.null(filematrix_colnames) && length(filematrix_colnames) == p) {
        rownames(rotation) <- filematrix_colnames
        names(mean_vec) <- filematrix_colnames
        names(column_sd) <- filematrix_colnames
        if (!is.null(scale_vec)) {
            names(scale_vec) <- filematrix_colnames
        }
        dimnames(covariance) <- list(filematrix_colnames, filematrix_colnames)
    }
    names(sdev) <- colnames(rotation)
    names(eigenvalues) <- colnames(rotation)
    names(explained) <- colnames(rotation)
    names(cumulative) <- colnames(rotation)

    scores <- NULL
    if (isTRUE(return_scores)) {
        scores <- pca_scores_stream_filematrix(
            X,
            rotation = rotation,
            center = center_vec,
            scale = scale_vec,
            ncomp = ncomp,
            chunk_size = provider$chunk_size
        )
    }

    result <- list(
        sdev = sdev,
        rotation = rotation,
        center = if (center) mean_vec else NULL,
        scale = scale_vec,
        scores = scores,
        column_sd = column_sd,
        eigenvalues = eigenvalues,
        explained_variance = explained,
        cumulative_variance = cumulative,
        covariance = covariance,
        nobs = n
    )
    new_bigpca_result(result, "stream_filematrix")
}

.spca_stream_row_provider <- function(provider,
                                      X,
                                      ncomp,
                                      center,
                                      scale,
                                      max_iter,
                                      tol,
                                      return_scores,
                                      verbose) {
    if (scale && !center) {
        stop("Scaling requires centring the variables", call. = FALSE)
    }

    n <- provider$nrow()
    p <- provider$ncol()
    if (n < 2L) {
        stop("PCA requires at least two observations", call. = FALSE)
    }

    max_iter <- .as_positive_scalar_integer(max_iter, "max_iter")
    tol <- as.numeric(tol)
    if (length(tol) != 1L || is.na(tol) || !is.finite(tol) || tol < 0) {
        stop("`tol` must be a non-negative finite number", call. = FALSE)
    }

    if (is.null(ncomp) || ncomp <= 0L) {
        ncomp <- min(n, p)
    }
    ncomp <- as.integer(min(ncomp, n, p))
    if (is.na(ncomp) || ncomp < 1L) {
        stop("`ncomp` must select at least one component", call. = FALSE)
    }

    block_starts <- seq.int(1L, n, by = provider$chunk_size)
    col_sum <- numeric(p)
    col_sum_sq <- numeric(p)
    observed <- 0L
    for (start in block_starts) {
        end <- min(start + provider$chunk_size - 1L, n)
        block <- provider$get_rows(start, end)
        .check_finite_numeric_block(block)
        col_sum <- col_sum + colSums(block)
        col_sum_sq <- col_sum_sq + colSums(block^2)
        observed <- observed + nrow(block)
    }
    if (observed != n) {
        stop("Provider returned an unexpected number of rows", call. = FALSE)
    }

    mean_vec <- col_sum / n
    col_ss <- numeric(p)
    for (start in block_starts) {
        end <- min(start + provider$chunk_size - 1L, n)
        block <- provider$get_rows(start, end)
        centered <- sweep(block, 2, mean_vec, "-", check.margin = FALSE)
        col_ss <- col_ss + colSums(centered^2)
    }
    column_sd <- sqrt(pmax(col_ss / max(1L, n - 1L), 0))
    center_vec <- if (center) mean_vec else NULL
    scale_vec <- if (scale) {
        out <- column_sd
        out[out == 0] <- 1
        out
    } else {
        NULL
    }

    total_variance <- if (scale) {
        as.numeric(p)
    } else if (center) {
        sum(column_sd^2)
    } else {
        sum(col_sum_sq) / n
    }

    transform_block <- function(block) {
        if (!is.null(center_vec)) {
            block <- sweep(block, 2, center_vec, "-", check.margin = FALSE)
        }
        if (!is.null(scale_vec)) {
            block <- sweep(block, 2, scale_vec, "/", check.margin = FALSE)
        }
        block
    }

    apply_operator <- function(Q_mat) {
        accum <- matrix(0, nrow = p, ncol = ncol(Q_mat))
        for (start in block_starts) {
            end <- min(start + provider$chunk_size - 1L, n)
            block <- transform_block(provider$get_rows(start, end))
            block_times_q <- block %*% Q_mat
            accum <- accum + crossprod(block, block_times_q)
        }
        accum
    }

    random_basis <- matrix(rnorm(p * ncomp), nrow = p, ncol = ncomp)
    Q <- svd(random_basis, nu = ncomp, nv = 0L)$u[, seq_len(ncomp), drop = FALSE]

    converged <- FALSE
    final_delta <- NA_real_
    iterations <- 0L
    for (iter in seq_len(max_iter)) {
        iterations <- iter
        operator_Q <- apply_operator(Q)
        Q_new <- svd(operator_Q, nu = ncomp, nv = 0L)$u[, seq_len(ncomp), drop = FALSE]
        final_delta <- base::norm(Q_new %*% t(Q_new) - Q %*% t(Q), type = "F")
        if (verbose) {
            message(sprintf("Iteration %d, subspace delta = %.6g", iter, final_delta))
        }
        Q <- Q_new
        if (final_delta <= tol) {
            converged <- TRUE
            break
        }
    }
    if (verbose && !converged) {
        message("Maximum iterations reached without meeting tolerance")
    }

    operator_Q <- apply_operator(Q)
    denom <- if (center) n - 1L else n
    if (denom <= 0L) {
        stop("Unable to compute streaming PCA with fewer than two observations", call. = FALSE)
    }
    small_matrix <- crossprod(Q, operator_Q) / denom
    eig <- eigen((small_matrix + t(small_matrix)) / 2, symmetric = TRUE)
    rotation <- Q %*% eig$vectors[, seq_len(ncomp), drop = FALSE]
    sdev <- sqrt(pmax(eig$values[seq_len(ncomp)], 0))
    eigenvalues <- sdev^2
    explained <- if (total_variance > 0) eigenvalues / total_variance else rep(0, length(eigenvalues))
    cumulative <- cumsum(explained)

    colnames(rotation) <- paste0("PC", seq_len(ncomp))
    filematrix_colnames <- .filematrix_colnames(X)
    if (!is.null(filematrix_colnames) && length(filematrix_colnames) == nrow(rotation)) {
        rownames(rotation) <- filematrix_colnames
        names(mean_vec) <- filematrix_colnames
        names(column_sd) <- filematrix_colnames
        if (!is.null(scale_vec)) {
            names(scale_vec) <- filematrix_colnames
        }
    }
    names(sdev) <- colnames(rotation)
    names(eigenvalues) <- colnames(rotation)
    names(explained) <- colnames(rotation)
    names(cumulative) <- colnames(rotation)

    scores <- NULL
    if (isTRUE(return_scores)) {
        scores <- pca_scores_stream_filematrix(
            X,
            rotation = rotation,
            center = center_vec,
            scale = scale_vec,
            ncomp = ncomp,
            chunk_size = provider$chunk_size
        )
    }

    result <- list(
        sdev = sdev,
        rotation = rotation,
        center = if (center) mean_vec else NULL,
        scale = scale_vec,
        scores = scores,
        column_sd = column_sd,
        eigenvalues = eigenvalues,
        explained_variance = explained,
        cumulative_variance = cumulative,
        covariance = NULL,
        nobs = n
    )
    result <- new_bigpca_result(result, "spca_filematrix")
    attr(result, "iterations") <- iterations
    attr(result, "tolerance") <- tol
    attr(result, "converged") <- isTRUE(converged)
    attr(result, "delta") <- final_delta
    result
}

.scores_stream_row_provider <- function(provider,
                                        X,
                                        rotation,
                                        center,
                                        scale,
                                        ncomp) {
    rotation <- as.matrix(rotation)
    if (!is.numeric(rotation)) {
        stop("`rotation` must be a numeric matrix", call. = FALSE)
    }
    if (nrow(rotation) != provider$ncol()) {
        stop("`rotation` must have one row per filematrix column", call. = FALSE)
    }
    if (ncol(rotation) < 1L) {
        stop("`rotation` must contain at least one component", call. = FALSE)
    }

    if (is.null(ncomp) || ncomp <= 0L) {
        ncomp <- ncol(rotation)
    }
    ncomp <- as.integer(ncomp)
    if (is.na(ncomp) || ncomp < 1L || ncomp > ncol(rotation)) {
        stop("`ncomp` must be between 1 and the number of rotation columns, or non-positive for all columns", call. = FALSE)
    }
    rotation <- rotation[, seq_len(ncomp), drop = FALSE]

    n <- provider$nrow()
    center_vec <- resolve_center_vector(center, provider$ncol())
    scale_vec <- resolve_scale_vector(scale, provider$ncol())
    block_starts <- seq.int(1L, n, by = provider$chunk_size)
    out <- matrix(0, nrow = n, ncol = ncomp)

    row_index <- 1L
    for (start in block_starts) {
        end <- min(start + provider$chunk_size - 1L, n)
        block <- provider$get_rows(start, end)
        .check_finite_numeric_block(block)
        block <- sweep(block, 2, center_vec, "-", check.margin = FALSE)
        block <- sweep(block, 2, scale_vec, "/", check.margin = FALSE)
        block_scores <- block %*% rotation
        rows <- seq.int(row_index, length.out = nrow(block))
        out[rows, ] <- block_scores
        row_index <- row_index + nrow(block)
    }

    score_colnames <- colnames(rotation)
    if (is.null(score_colnames)) {
        score_colnames <- paste0("PC", seq_len(ncomp))
    }
    row_names <- .filematrix_rownames(X)

    colnames(out) <- score_colnames
    if (!is.null(row_names) && length(row_names) == n) {
        rownames(out) <- row_names
    }
    out
}

.check_finite_numeric_block <- function(block) {
    if (!is.numeric(block)) {
        stop("Provider block is not numeric; unsupported filematrix storage.", call. = FALSE)
    }
    if (any(!is.finite(block))) {
        stop("Provider block contains non-finite values.", call. = FALSE)
    }
    invisible(TRUE)
}

.check_filematrix_covariance_budget <- function(p) {
    max_gb <- getOption("bigPCAcpp.filematrix.max_cov_gb", 2)
    if (is.null(max_gb) || identical(max_gb, Inf)) {
        return(invisible(TRUE))
    }
    max_gb <- suppressWarnings(as.numeric(max_gb))
    if (length(max_gb) != 1L || is.na(max_gb) || max_gb < 0) {
        stop("Option `bigPCAcpp.filematrix.max_cov_gb` must be a non-negative numeric scalar or Inf.", call. = FALSE)
    }
    cov_gb <- (as.numeric(p) * as.numeric(p) * 8) / 1024^3
    if (cov_gb > max_gb) {
        stop(
            sprintf(
                "Refusing exact filematrix PCA: the %d x %d covariance matrix would require %.3f GB, exceeding option `bigPCAcpp.filematrix.max_cov_gb` = %.3f GB. Use pca_spca_stream_filematrix() for very-wide filematrix workflows or increase the option for benchmark validation.",
                p,
                p,
                cov_gb,
                max_gb
            ),
            call. = FALSE
        )
    }
    invisible(TRUE)
}

.filematrix_colnames <- function(X) {
    tryCatch(colnames(X), error = function(e) NULL)
}

.filematrix_rownames <- function(X) {
    tryCatch(rownames(X), error = function(e) NULL)
}
