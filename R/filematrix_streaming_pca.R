#' Experimental filematrix streaming PCA helpers
#'
#' @description
#' Experimental streaming helpers for `filematrix` inputs. These functions read
#' the source matrix in row blocks through [make_filematrix_row_provider()],
#' avoiding `bigmemory` memory-mapped access patterns on shared filesystems.
#' They are intended for streaming SPCA and score projection workflows, not as
#' replacements for exact covariance PCA or highly optimised random-access
#' workloads.
#'
#' @inheritParams pca_spca
#' @param X A `filematrix` object.
#' @param chunk_size Number of rows to read per block.
#' @param sink For `pca_scores_stream_filematrix()`, where scores should be
#'   written. `"r"` returns a dense R matrix, `"none"` scans and validates
#'   without retaining scores, and `"filematrix"` writes to a new `filematrix`.
#' @param filematrix_base Optional filename base used when `sink = "filematrix"`.
#' @param type Storage type for the output scores when `sink = "filematrix"`.
#'
#' @return `pca_spca_stream_filematrix()` returns a [`bigpca`] object compatible
#'   with [pca_spca()] outputs. `pca_scores_stream_filematrix()` returns a dense
#'   matrix for `sink = "r"`, a `filematrix` object for `sink = "filematrix"`,
#'   and an invisible scan summary for `sink = "none"`.
#'
#' @name filematrix_streaming_pca
NULL

#' @rdname filematrix_streaming_pca
#' @export
pca_spca_stream_filematrix <- function(
        X,
        ncomp = NULL,
        center = TRUE,
        scale = FALSE,
        chunk_size = 8192L,
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
#' @param center Optional numeric vector of column centers.
#' @param scale Optional numeric vector of column scale factors.
#' @export
pca_scores_stream_filematrix <- function(
        X,
        rotation,
        center = NULL,
        scale = NULL,
        chunk_size = 8192L,
        sink = c("r", "filematrix", "none"),
        filematrix_base = NULL,
        type = "double") {
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
        sink = sink,
        filematrix_base = filematrix_base,
        type = type
    )
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
            chunk_size = provider$chunk_size,
            sink = "r"
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
    result <- new_bigpca_result(result, "spca_stream_filematrix")
    attr(result, "iterations") <- iterations
    attr(result, "tolerance") <- tol
    attr(result, "converged") <- isTRUE(converged)
    attr(result, "delta") <- final_delta
    attr(result, "storage_type") <- provider$storage_type
    attr(result, "experimental") <- TRUE
    result
}

.scores_stream_row_provider <- function(provider,
                                        X,
                                        rotation,
                                        center,
                                        scale,
                                        sink,
                                        filematrix_base,
                                        type) {
    sink <- match.arg(sink, c("r", "filematrix", "none"))
    rotation <- as.matrix(rotation)
    if (!is.numeric(rotation)) {
        stop("`rotation` must be a numeric matrix", call. = FALSE)
    }
    if (nrow(rotation) != provider$ncol()) {
        stop("`rotation` must have one row per filematrix column", call. = FALSE)
    }

    n <- provider$nrow()
    k <- ncol(rotation)
    center_vec <- resolve_center_vector(center, provider$ncol())
    scale_vec <- resolve_scale_vector(scale, provider$ncol())
    block_starts <- seq.int(1L, n, by = provider$chunk_size)

    out <- switch(
        sink,
        r = matrix(0, nrow = n, ncol = k),
        filematrix = {
            if (is.null(filematrix_base)) {
                filematrix_base <- tempfile("bigPCAcpp_filematrix_scores_")
            }
            filematrix::fm.create(
                filenamebase = filematrix_base,
                nrow = n,
                ncol = k,
                type = type
            )
        },
        none = NULL
    )

    row_index <- 1L
    for (start in block_starts) {
        end <- min(start + provider$chunk_size - 1L, n)
        block <- provider$get_rows(start, end)
        .check_finite_numeric_block(block)
        block <- sweep(block, 2, center_vec, "-", check.margin = FALSE)
        block <- sweep(block, 2, scale_vec, "/", check.margin = FALSE)
        block_scores <- block %*% rotation
        rows <- seq.int(row_index, length.out = nrow(block))
        if (sink == "r") {
            out[rows, ] <- block_scores
        } else if (sink == "filematrix") {
            out[rows, ] <- block_scores
        }
        row_index <- row_index + nrow(block)
    }

    score_colnames <- colnames(rotation)
    if (is.null(score_colnames)) {
        score_colnames <- paste0("PC", seq_len(k))
    }
    row_names <- .filematrix_rownames(X)

    if (sink == "none") {
        return(invisible(list(
            nobs = n,
            ncomp = k,
            sink = "none",
            storage_type = provider$storage_type
        )))
    }

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

.filematrix_colnames <- function(X) {
    tryCatch(colnames(X), error = function(e) NULL)
}

.filematrix_rownames <- function(X) {
    tryCatch(rownames(X), error = function(e) NULL)
}
