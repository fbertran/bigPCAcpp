#' Robust principal component analysis
#'
#' @description
#' Compute principal component analysis (PCA) using robust measures of location
#' and scale so that extreme observations have a reduced influence on the
#' resulting components. The implementation centres each variable by its median
#' and, when requested, scales by the median absolute deviation (MAD) before
#' performing an iteratively reweighted singular value decomposition that
#' down-weights observations with unusually large reconstruction errors.
#'
#' @param x A numeric matrix, data frame, or an object coercible to a numeric
#'   matrix. Missing values are not supported.
#' @param center Logical; should variables be centred by their median before
#'   applying PCA?
#' @param scale Logical; when `TRUE`, variables are scaled by the MAD after
#'   centring. Scaling requires `center = TRUE`.
#' @param ncomp Number of components to retain. Use `NULL` or a non-positive
#'   value to keep all components returned by the decomposition.
#'
#' @return A [`bigpca`] object mirroring the structure of
#'   [pca_bigmatrix()] with robust estimates of location, scale, and
#'   variance metrics.
#'
#' @examples
#' set.seed(42)
#' x <- matrix(rnorm(50), nrow = 10)
#' x[1, 1] <- 25  # outlier
#' robust <- pca_robust(x, ncomp = 2)
#' robust$sdev
#'
#' @export
pca_robust <- function(x, center = TRUE, scale = FALSE, ncomp = NULL) {
    if (!is.matrix(x)) {
        if (methods::is(x, "big.matrix")) {
            if (!requireNamespace("bigmemory", quietly = TRUE)) {
                stop("`bigmemory` must be installed to coerce big.matrix inputs", call. = FALSE)
            }
            x <- bigmemory::as.matrix(x)
        } else {
            x <- as.matrix(x)
        }
    }
    if (!is.numeric(x)) {
        stop("`x` must contain numeric values for PCA", call. = FALSE)
    }
    if (anyNA(x)) {
        stop("`x` must not contain missing values for robust PCA", call. = FALSE)
    }
    n <- nrow(x)
    p <- ncol(x)
    if (n == 0 || p == 0) {
        stop("`x` must contain at least one observation and one variable", call. = FALSE)
    }
    if (!center && scale) {
        stop("Robust scaling requires `center = TRUE`", call. = FALSE)
    }
    if (is.null(ncomp) || ncomp <= 0) {
        ncomp <- min(n, p)
    } else {
        ncomp <- min(ncomp, n, p)
    }

    center_vec <- if (center) {
        apply(x, 2, stats::median)
    } else {
        rep.int(0, p)
    }
    centered <- sweep(x, 2, center_vec, "-", check.margin = FALSE)

    if (scale) {
        scale_vec <- apply(centered, 2, stats::mad)
        scale_vec[!is.finite(scale_vec) | scale_vec == 0] <- 1
        scaled <- sweep(centered, 2, scale_vec, "/", check.margin = FALSE)
    } else {
        scale_vec <- NULL
        scaled <- centered
    }

    svd_res <- svd_robust(scaled, ncomp = ncomp)
    kept <- seq_len(ncomp)
    rotation <- svd_res$v[, kept, drop = FALSE]
    base_weights <- as.numeric(svd_res$weights)
    weight_scale <- if (length(base_weights) == n) {
        sqrt(pmax(base_weights, 1e-6))
    } else {
        rep(1, n)
    }
    weighted_scores <- sweep(svd_res$u[, kept, drop = FALSE], 2, svd_res$d[kept], `*`)
    scores <- sweep(weighted_scores, 1, weight_scale, "/")
    leverage_weights <- {
        score_norm <- sqrt(rowSums(scores^2))
        scale_norm <- stats::mad(score_norm, center = 0, constant = 1)
        if (!is.finite(scale_norm) || scale_norm <= .Machine$double.eps) {
            scale_norm <- sqrt(mean(score_norm^2))
        }
        if (!is.finite(scale_norm) || scale_norm <= .Machine$double.eps) {
            rep(1, length(score_norm))
        } else {
            scaled_norm <- score_norm / (scale_norm * 1.345)
            weights <- ifelse(!is.finite(scaled_norm) | scaled_norm <= 1, 1, 1 / scaled_norm)
            pmin(pmax(weights, 1e-6), 1)
        }
    }
    if (length(base_weights) && length(base_weights) == length(leverage_weights)) {
        combined_weights <- pmin(base_weights, leverage_weights)
    } else {
        combined_weights <- leverage_weights
    }
    weight_contrib <- if (length(combined_weights) == n) {
        combined_weights^2
    } else {
        rep(1, n)
    }
    weighted_norms <- colSums(sweep(scores^2, 1, weight_contrib, `*`))
    total_weight <- sum(weight_contrib)
    singular_values <- sqrt(weighted_norms)
    sdev <- if (total_weight > 1) {
        sqrt(weighted_norms / (total_weight - 1))
    } else {
        singular_values
    }

    column_sd <- apply(centered, 2, stats::mad)
    eigenvalues <- sdev^2
    total_var <- sum(eigenvalues)
    explained <- if (total_var > 0) eigenvalues / total_var else rep(0, length(eigenvalues))
    cumulative <- cumsum(explained)

    result <- list(
        sdev = sdev,
        rotation = rotation,
        center = if (center) center_vec else NULL,
        scale = if (scale) scale_vec else NULL,
        scores = scores,
        column_sd = column_sd,
        eigenvalues = eigenvalues,
        explained_variance = explained,
        cumulative_variance = cumulative,
        nobs = n,
        robust_weights = combined_weights,
        robust_iterations = svd_res$iterations
    )
    new_bigpca_result(result, "robust")
}

#' Prepare iteratively reweighted singular value decomposition
#'
#' @description
#' Internal helper used by [pca_robust()] to compute a singular value
#' decomposition that is less sensitive to individual rows with extreme values.
#' The routine alternates between computing the SVD of a row-weighted matrix and
#' updating the weights via a Huber-type scheme based on the reconstruction
#' residuals.
#'
#' @param x Numeric matrix for which the decomposition should be computed.
#' @param ncomp Number of leading components to retain.
#' @param max_iter Maximum number of reweighting iterations.
#' @param tol Convergence tolerance applied to successive changes in the row
#'   weights and singular values.
#' @param huber_k Tuning constant controlling the aggressiveness of the Huber
#'   weight function. Larger values down-weight fewer observations.
#'
#' @return A list containing x, n, p, ncomp, max_iter, tol and huber_k.
#' @keywords internal
prepare_svd_robust_input <- function(x, ncomp, max_iter, tol, huber_k) {
    if (!is.matrix(x)) {
        stop("`x` must be a matrix for robust SVD", call. = FALSE)
    }
    if (!is.numeric(x)) {
        stop("`x` must contain numeric values for robust SVD", call. = FALSE)
    }
    n <- nrow(x)
    p <- ncol(x)
    if (n == 0L || p == 0L) {
        stop("`x` must contain at least one row and one column", call. = FALSE)
    }
    if (is.null(ncomp) || ncomp <= 0) {
        ncomp <- min(n, p)
    } else {
        ncomp <- min(ncomp, n, p)
    }
    max_iter <- as.integer(max_iter)
    if (length(max_iter) != 1L || is.na(max_iter) || max_iter <= 0L) {
        stop("`max_iter` must be a positive integer", call. = FALSE)
    }
    tol <- as.numeric(tol)
    if (length(tol) != 1L || !is.finite(tol) || tol <= 0) {
        stop("`tol` must be a positive finite number", call. = FALSE)
    }
    huber_k <- as.numeric(huber_k)
    if (length(huber_k) != 1L || !is.finite(huber_k) || huber_k <= 0) {
        stop("`huber_k` must be a positive finite number", call. = FALSE)
    }
    list(
        x = x,
        n = n,
        p = p,
        ncomp = ncomp,
        max_iter = max_iter,
        tol = tol,
        huber_k = huber_k
    )
}

#' Iteratively reweighted singular value decomposition
#'
#' @description Internal helper used by [pca_robust()] to compute a singular value
#' decomposition that is less sensitive to individual rows with extreme values.
#' The routine alternates between computing the SVD of a row-weighted matrix and
#' updating the weights via a Huber-type scheme based on the reconstruction
#' residuals.
#'
#' @inheritParams svd_robust
#' @return A list containing the left and right singular vectors (`u` and `v`),
#'   the singular values (`d`), the final row weights (`weights`), and the number
#'   of iterations required for convergence (`iterations`). The structure mirrors
#'   base R's [base::svd()] output with additional metadata.
#'   
svd_robust_R <- function(x, ncomp, max_iter = 25L, tol = sqrt(.Machine$double.eps), huber_k = 1.345) {
  prepared <- prepare_svd_robust_input(x, ncomp, max_iter, tol, huber_k)
  x <- prepared$x
  ncomp <- prepared$ncomp
  max_iter <- prepared$max_iter
  tol <- prepared$tol
  huber_k <- prepared$huber_k
  
  n <- nrow(x)
  p <- ncol(x)
  nu <- min(n, ncomp)
  nv <- min(p, ncomp)
  weights <- rep(1, n)
  sqrt_weights <- rep(1, n)
  prev_singular <- rep(0, ncomp)
  last_svd <- base::svd(sweep(x, 1, sqrt_weights, `*`), nu = nu, nv = nv)
  iterations <- 1L
  
  for (iter in seq_len(max_iter)) {
    if (iter > 1L) {
      weighted <- sweep(x, 1, sqrt_weights, `*`)
      last_svd <- base::svd(weighted, nu = nu, nv = nv)
      iterations <- iter
    }
    
    kept <- seq_len(ncomp)
    rotation <- last_svd$v[, kept, drop = FALSE]
    u_trunc <- last_svd$u[, kept, drop = FALSE]
    singular <- last_svd$d[kept]
    
    scores_weighted <- sweep(u_trunc, 2, singular, `*`)
    fitted_weighted <- scores_weighted %*% t(rotation)
    fitted <- sweep(fitted_weighted, 1, sqrt_weights, `/`)
    residual <- x - fitted
    row_resid <- sqrt(rowSums(residual^2))
    
    scale <- stats::mad(row_resid, center = 0, constant = 1)
    if (!is.finite(scale) || scale <= .Machine$double.eps) {
      scale <- sqrt(mean(row_resid^2))
    }
    if (!is.finite(scale) || scale <= .Machine$double.eps) {
      break
    }
    
    scaled <- row_resid / (scale * huber_k)
    new_weights <- ifelse(!is.finite(scaled) | scaled <= 1, 1, 1 / scaled)
    new_weights[!is.finite(new_weights)] <- 1
    new_weights <- pmin(pmax(new_weights, 1e-6), 1)
    
    weight_change <- max(abs(new_weights - weights))
    rel_change <- if (iter == 1L) Inf else {
      denom <- max(1, max(prev_singular))
      max(abs(singular - prev_singular)) / denom
    }
    
    weights <- new_weights
    sqrt_weights <- sqrt(weights)
    prev_singular <- singular
    
    if (!is.finite(weight_change)) {
      weight_change <- 0
    }
    if (weight_change < tol && rel_change < tol) {
      iterations <- iter
      break
    }
    iterations <- iter
  }
  
  list(
    u = last_svd$u[, seq_len(nu), drop = FALSE],
    d = last_svd$d[seq_len(ncomp)],
    v = last_svd$v[, seq_len(nv), drop = FALSE],
    weights = weights,
    iterations = iterations
  )
}


#' Robust singular value decomposition (C++ backend)
#'
#' @description
#' Compute the iteratively reweighted SVD using the high-performance C++
#' implementation. The interface mirrors [svd_robust_R()] while delegating the
#' heavy lifting to compiled code.
#'
#' @param x Numeric matrix for which the decomposition should be computed.
#' @param ncomp Number of leading components to retain.
#' @param max_iter Maximum number of reweighting iterations.
#' @param tol Convergence tolerance applied to successive changes in the row
#'   weights and singular values.
#' @param huber_k Tuning constant controlling the aggressiveness of the Huber
#'   weight function. Larger values down-weight fewer observations.
#' @return A list containing the left and right singular vectors (`u` and `v`),
#'   the singular values (`d`), the final row weights (`weights`), and the number
#'   of iterations required for convergence (`iterations`).
#' @keywords internal
svd_robust <- function(x, ncomp, max_iter = 25L, tol = sqrt(.Machine$double.eps), huber_k = 1.345) {
    prepared <- prepare_svd_robust_input(x, ncomp, max_iter, tol, huber_k)
    .svd_robust_cpp(prepared$x, prepared$ncomp, prepared$max_iter, prepared$tol, prepared$huber_k)
}




#' Scalable principal component analysis via streaming power iterations
#'
#' @description
#' Implements the scalable PCA (sPCA) procedure of Elgamal et al. (2015), which
#' uses block power iterations to approximate the leading principal components
#' while streaming the data in manageable chunks. The algorithm only requires
#' matrix-vector products, allowing large matrices to be processed without
#' materialising the full cross-product in memory.
#'
#' @param x A numeric matrix, data frame, or [`bigmemory::big.matrix`]. The
#'   input is processed in row-wise blocks so that large matrices can be
#'   analysed without creating dense copies in R memory.
#' @param ncomp Number of principal components to retain. Use `NULL` or a
#'   non-positive value to keep `min(nrow(x), ncol(x))` components.
#' @param center Logical; should column means be subtracted before performing
#'   PCA?
#' @param scale Logical; when `TRUE`, columns are scaled to unit variance after
#'   centring. Scaling requires `center = TRUE`.
#' @param block_size Number of rows to stream per block when computing column
#'   statistics and matrix-vector products.
#' @param max_iter Maximum number of block power iterations.
#' @param tol Convergence tolerance applied to the Frobenius norm of the
#'   difference between successive subspace projectors.
#' @param seed Optional integer seed used to initialise the random starting
#'   basis.
#' @param return_scores Logical; when `TRUE`, principal component scores are
#'   computed in a final streaming pass over the data.
#' @param verbose Logical; when `TRUE`, diagnostic messages describing the
#'   iteration progress are emitted.
#'
#' @return A [`bigpca`] object containing the approximate PCA solution with the
#'   same structure as [pca_bigmatrix()]. The result includes component standard
#'   deviations, rotation/loadings, optional scores, column statistics, and
#'   variance summaries. Additional metadata is stored in
#'   `attr(result, "iterations")` (number of iterations performed),
#'   `attr(result, "tolerance")` (requested tolerance), and
#'   `attr(result, "converged")` (logical convergence flag).
#'
#' @references
#' Tarek Elgamal, Maysam Yabandeh, Ashraf Aboulnaga, Waleed Mustafa, and Mohamed
#' Hefeeda (2015). *sPCA: Scalable Principal Component Analysis for Big Data on
#' Distributed Platforms*. Proceedings of the 2015 ACM SIGMOD International
#' Conference on Management of Data. <doi:10.1145/2723372.2751520>.
#'
#' @export
pca_spca <- function(x,
                     ncomp = NULL,
                     center = TRUE,
                     scale = FALSE,
                     block_size = 2048L,
                     max_iter = 50L,
                     tol = 1e-4,
                     seed = NULL,
                     return_scores = FALSE,
                     verbose = FALSE) {
    if (!is.null(seed)) {
        old_seed <- if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)) {
            get(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
        } else {
            NULL
        }
        on.exit({
            if (is.null(old_seed)) {
                if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)) {
                    rm(".Random.seed", envir = .GlobalEnv)
                }
            } else {
                assign(".Random.seed", old_seed, envir = .GlobalEnv)
            }
        }, add = TRUE)
        set.seed(seed)
    }

    if (scale && !center) {
        stop("Scaling requires centring the variables", call. = FALSE)
    }

    if (inherits(x, "big.matrix")) {
        if (!requireNamespace("bigmemory", quietly = TRUE)) {
            stop("`bigmemory` must be installed to work with big.matrix inputs", call. = FALSE)
        }
        accessor <- list(
            n = nrow(x),
            p = ncol(x),
            colnames = colnames(x),
            rownames = rownames(x),
            pull = function(start, end) {
                bigmemory::as.matrix(x[start:end, , drop = FALSE])
            }
        )
    } else {
        if (!is.matrix(x)) {
            x <- as.matrix(x)
        }
        if (!is.numeric(x)) {
            stop("`x` must contain numeric values for scalable PCA", call. = FALSE)
        }
        accessor <- list(
            n = nrow(x),
            p = ncol(x),
            colnames = colnames(x),
            rownames = rownames(x),
            pull = function(start, end) {
                x[start:end, , drop = FALSE]
            }
        )
    }

    n <- accessor$n
    p <- accessor$p
    if (n == 0 || p == 0) {
        stop("`x` must contain at least one observation and one variable", call. = FALSE)
    }
    if (n < 2) {
        stop("PCA requires at least two observations", call. = FALSE)
    }

    block_size <- as.integer(block_size)
    if (is.na(block_size) || block_size <= 0L) {
        stop("`block_size` must be a positive integer", call. = FALSE)
    }
    block_size <- min(block_size, n)

    max_iter <- as.integer(max_iter)
    if (is.na(max_iter) || max_iter <= 0L) {
        stop("`max_iter` must be a positive integer", call. = FALSE)
    }

    if (is.null(ncomp) || ncomp <= 0L) {
        ncomp <- min(n, p)
    }
    ncomp <- as.integer(min(ncomp, n, p))
    if (ncomp < 1L) {
        stop("`ncomp` must select at least one component", call. = FALSE)
    }

    block_starts <- seq.int(1L, n, by = block_size)
    col_sum <- numeric(p)
    col_sum_sq <- numeric(p)
    observed <- 0L
    for (start in block_starts) {
        end <- min(start + block_size - 1L, n)
        block <- accessor$pull(start, end)
        if (!is.numeric(block)) {
            stop("`x` must contain numeric values for scalable PCA", call. = FALSE)
        }
        col_sum <- col_sum + colSums(block)
        col_sum_sq <- col_sum_sq + colSums(block^2)
        observed <- observed + nrow(block)
    }
    if (observed != n) {
        n <- observed
        block_size <- min(block_size, n)
        block_starts <- seq.int(1L, n, by = block_size)
    }

    mean_vec <- col_sum / n
    col_ss <- numeric(p)
    for (start in block_starts) {
        end <- min(start + block_size - 1L, n)
        block <- accessor$pull(start, end)
        centered <- sweep(block, 2, mean_vec, "-", check.margin = FALSE)
        col_ss <- col_ss + colSums(centered^2)
    }
    denom_sd <- max(1, n - 1L)
    column_sd <- sqrt(pmax(col_ss / denom_sd, 0))
    scale_vec <- if (scale) {
        out <- column_sd
        out[out == 0] <- 1
        out
    } else {
        NULL
    }

    center_vec <- if (center) mean_vec else NULL
    total_variance <- if (scale) {
        if (center) {
            as.numeric(p)
        } else {
            scale_sq <- if (!is.null(scale_vec)) scale_vec^2 else rep(1, length(col_sum_sq))
            sum(col_sum_sq / scale_sq) / n
        }
    } else {
        if (center) {
            sum(column_sd^2)
        } else {
            sum(col_sum_sq) / n
        }
    }

    random_basis <- matrix(rnorm(p * ncomp), nrow = p, ncol = ncomp)
    svd_init <- svd(random_basis, nu = ncomp, nv = 0L)
    Q <- svd_init$u[, seq_len(ncomp), drop = FALSE]

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
            end <- min(start + block_size - 1L, n)
            block <- transform_block(accessor$pull(start, end))
            if (nrow(block) == 0L) {
                next
            }
            block_times_q <- block %*% Q_mat
            accum <- accum + crossprod(block, block_times_q)
        }
        accum
    }

    convergence_flag <- FALSE
    final_delta <- NA_real_
    iterations <- 0L
    for (iter in seq_len(max_iter)) {
        iterations <- iter
        operator_Q <- apply_operator(Q)
        svd_res <- svd(operator_Q, nu = ncomp, nv = 0L)
        Q_new <- svd_res$u[, seq_len(ncomp), drop = FALSE]
        projector_diff <- Q_new %*% t(Q_new) - Q %*% t(Q)
        final_delta <- base::norm(projector_diff, type = "F")
        if (verbose) {
            message(sprintf("Iteration %d, subspace delta = %.6g", iter, final_delta))
        }
        Q <- Q_new
        if (final_delta <= tol) {
            convergence_flag <- TRUE
            break
        }
    }
    if (verbose && !convergence_flag) {
        message("Maximum iterations reached without meeting tolerance")
    }

    operator_Q <- apply_operator(Q)
    denom <- if (center) n - 1 else n
    if (denom <= 0) {
        stop("Unable to compute covariance with fewer than two observations", call. = FALSE)
    }
    small_matrix <- crossprod(Q, operator_Q) / denom
    sym_small <- (small_matrix + t(small_matrix)) / 2
    eig <- eigen(sym_small, symmetric = TRUE)
    rotation <- Q %*% eig$vectors[, seq_len(ncomp), drop = FALSE]
    sdev <- sqrt(pmax(eig$values[seq_len(ncomp)], 0))
    eigenvalues <- sdev^2

    explained <- if (total_variance > 0) eigenvalues / total_variance else rep(0, length(eigenvalues))
    cumulative <- cumsum(explained)

    if (!is.null(accessor$colnames)) {
        rownames(rotation) <- accessor$colnames
    }
    colnames(rotation) <- paste0("PC", seq_len(ncomp))
    names(sdev) <- colnames(rotation)
    names(eigenvalues) <- colnames(rotation)
    names(explained) <- colnames(rotation)
    names(cumulative) <- colnames(rotation)

    scores <- NULL
    if (isTRUE(return_scores)) {
        scores <- matrix(0, nrow = n, ncol = ncomp)
        row_index <- 1L
        for (start in block_starts) {
            end <- min(start + block_size - 1L, n)
            block <- transform_block(accessor$pull(start, end))
            if (nrow(block) == 0L) {
                next
            }
            block_scores <- block %*% rotation
            rows <- seq.int(row_index, length.out = nrow(block))
            scores[rows, ] <- block_scores
            row_index <- row_index + nrow(block)
        }
        colnames(scores) <- colnames(rotation)
        if (!is.null(accessor$rownames)) {
            rownames(scores) <- accessor$rownames
        }
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
    result <- new_bigpca_result(result, "spca")
    attr(result, "iterations") <- iterations
    attr(result, "tolerance") <- tol
    attr(result, "converged") <- isTRUE(convergence_flag)
    attr(result, "delta") <- final_delta
    result
}

#' Principal component analysis for `bigmemory::big.matrix` inputs
#' 
#' @description
#' Implements the scalable PCA (sPCA) procedure of Elgamal et al. (2015), which
#' uses block power iterations to approximate the leading principal components
#' while streaming the data in manageable chunks. The algorithm only requires
#' matrix-vector products, allowing large matrices to be processed without
#' materialising the full cross-product in memory.
#'
#' @param x A numeric matrix, data frame, [`bigmemory::big.matrix`], or an
#'   external pointer referencing a big.matrix. The input is processed in
#'   row-wise blocks so that large matrices can be analysed without creating
#'   dense copies in R memory.
#' @param ncomp Number of principal components to retain. Use `NULL` or a
#'   non-positive value to keep `min(nrow(x), ncol(x))` components.
#' @param center Logical; should column means be subtracted before performing
#'   PCA?
#' @param scale Logical; when `TRUE`, columns are scaled to unit variance after
#'   centring. Scaling requires `center = TRUE`.
#' @param block_size Number of rows to stream per block when computing column
#'   statistics and matrix-vector products.
#' @param max_iter Maximum number of block power iterations.
#' @param tol Convergence tolerance applied to the Frobenius norm of the
#'   difference between successive subspace projectors.
#' @param seed Optional integer seed used to initialise the random starting
#'   basis.
#' @param return_scores Logical; when `TRUE`, principal component scores are
#'   computed in a final streaming pass over the data.
#' @param verbose Logical; when `TRUE`, diagnostic messages describing the
#'   iteration progress are emitted.
#'
#' @return A [`bigpca`] object containing the approximate PCA solution with the
#'   same structure as [pca_bigmatrix()]. The result includes component standard
#'   deviations, rotation/loadings, optional scores, column statistics, and
#'   variance summaries. Additional metadata is stored in
#'   `attr(result, "iterations")` (number of iterations performed),
#'   `attr(result, "tolerance")` (requested tolerance), and
#'   `attr(result, "converged")` (logical convergence flag).
#'
#' @references
#' Tarek Elgamal, Maysam Yabandeh, Ashraf Aboulnaga, Waleed Mustafa, and Mohamed
#' Hefeeda (2015). *sPCA: Scalable Principal Component Analysis for Big Data on
#' Distributed Platforms*. Proceedings of the 2015 ACM SIGMOD International
#' Conference on Management of Data. <doi:10.1145/2723372.2751520>.
#'
#' @name pca_spca
#' @export
pca_spca <- function(x,
                     ncomp = NULL,
                     center = TRUE,
                     scale = FALSE,
                     block_size = 2048L,
                     max_iter = 50L,
                     tol = 1e-4,
                     seed = NULL,
                     return_scores = FALSE,
                     verbose = FALSE) {
    if (scale && !center) {
        stop("Scaling requires centring the variables", call. = FALSE)
    }

    run_spca <- function() {
        if (inherits(x, "big.matrix") || typeof(x) == "externalptr") {
            ptr <- resolve_big_pointer(x, "x")
            block_size_int <- as.integer(block_size)
            if (is.na(block_size_int) || block_size_int <= 0L) {
                stop("`block_size` must be a positive integer", call. = FALSE)
            }
            max_iter_int <- as.integer(max_iter)
            if (is.na(max_iter_int) || max_iter_int <= 0L) {
                stop("`max_iter` must be a positive integer", call. = FALSE)
            }
            ncomp_int <- if (is.null(ncomp) || ncomp <= 0L) {
                -1L
            } else {
                as.integer(ncomp)
            }
            if (is.na(ncomp_int)) {
                stop("`ncomp` must be numeric", call. = FALSE)
            }
            result <- .pca_spca_bigmatrix(
                ptr,
                center,
                scale,
                ncomp_int,
                block_size_int,
                max_iter_int,
                tol,
                return_scores,
                verbose
            )
            result <- new_bigpca_result(result, "spca_bigmemory")
            rotation <- result$rotation
            comps <- seq_len(ncol(rotation))
            colnames(rotation) <- paste0("PC", comps)
            if (inherits(x, "big.matrix")) {
                cn <- tryCatch(colnames(x), error = function(e) NULL)
                rn <- tryCatch(rownames(x), error = function(e) NULL)
                if (!is.null(cn) && length(cn) == nrow(rotation)) {
                    rownames(rotation) <- cn
                }
                if (isTRUE(return_scores) && !is.null(result$scores) &&
                    !is.null(rn) && length(rn) == nrow(result$scores)) {
                    rownames(result$scores) <- rn
                }
            }
            result$rotation <- rotation
            if (!is.null(result$sdev)) {
                names(result$sdev) <- colnames(rotation)
            }
            if (!is.null(result$eigenvalues)) {
                names(result$eigenvalues) <- colnames(rotation)
            }
            if (!is.null(result$explained_variance)) {
                names(result$explained_variance) <- colnames(rotation)
            }
            if (!is.null(result$cumulative_variance)) {
                names(result$cumulative_variance) <- colnames(rotation)
            }
            if (isTRUE(return_scores) && !is.null(result$scores)) {
                colnames(result$scores) <- colnames(rotation)
            }
            return(result)
        }
        pca_spca_R(
            x = x,
            ncomp = ncomp,
            center = center,
            scale = scale,
            block_size = block_size,
            max_iter = max_iter,
            tol = tol,
            seed = NULL,
            return_scores = return_scores,
            verbose = verbose
        )
    }

    if (is.null(seed)) {
        run_spca()
    } else {
        old_seed <- if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)) {
            get(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
        } else {
            NULL
        }
        on.exit({
            if (is.null(old_seed)) {
                if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)) {
                    rm(".Random.seed", envir = .GlobalEnv)
                }
            } else {
                assign(".Random.seed", old_seed, envir = .GlobalEnv)
            }
        }, add = TRUE)
        set.seed(seed)
        run_spca()
    }
}

#' @rdname pca_spca
#' @export
pca_spca_R <- function(x,
                       ncomp = NULL,
                       center = TRUE,
                       scale = FALSE,
                       block_size = 2048L,
                       max_iter = 50L,
                       tol = 1e-4,
                       seed = NULL,
                       return_scores = FALSE,
                       verbose = FALSE) {
    if (!is.null(seed)) {
        old_seed <- if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)) {
            get(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
        } else {
            NULL
        }
        on.exit({
            if (is.null(old_seed)) {
                if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)) {
                    rm(".Random.seed", envir = .GlobalEnv)
                }
            } else {
                assign(".Random.seed", old_seed, envir = .GlobalEnv)
            }
        }, add = TRUE)
        set.seed(seed)
    }

    if (scale && !center) {
        stop("Scaling requires centring the variables", call. = FALSE)
    }

    if (!is.matrix(x)) {
        x <- as.matrix(x)
    }
    if (!is.numeric(x)) {
        stop("`x` must contain numeric values for scalable PCA", call. = FALSE)
    }

    n <- nrow(x)
    p <- ncol(x)
    if (n == 0 || p == 0) {
        stop("`x` must contain at least one observation and one variable", call. = FALSE)
    }
    if (n < 2) {
        stop("PCA requires at least two observations", call. = FALSE)
    }

    block_size <- as.integer(block_size)
    if (is.na(block_size) || block_size <= 0L) {
        stop("`block_size` must be a positive integer", call. = FALSE)
    }
    block_size <- min(block_size, n)

    max_iter <- as.integer(max_iter)
    if (is.na(max_iter) || max_iter <= 0L) {
        stop("`max_iter` must be a positive integer", call. = FALSE)
    }

    if (is.null(ncomp) || ncomp <= 0L) {
        ncomp <- min(n, p)
    }
    ncomp <- as.integer(min(ncomp, n, p))
    if (ncomp < 1L) {
        stop("`ncomp` must select at least one component", call. = FALSE)
    }

    block_starts <- seq.int(1L, n, by = block_size)
    col_sum <- numeric(p)
    col_sum_sq <- numeric(p)
    observed <- 0L
    for (start in block_starts) {
        end <- min(start + block_size - 1L, n)
        block <- x[start:end, , drop = FALSE]
        col_sum <- col_sum + colSums(block)
        col_sum_sq <- col_sum_sq + colSums(block^2)
        observed <- observed + nrow(block)
    }
    if (observed != n) {
        n <- observed
        block_size <- min(block_size, n)
        block_starts <- seq.int(1L, n, by = block_size)
    }

    mean_vec <- col_sum / n
    col_ss <- numeric(p)
    for (start in block_starts) {
        end <- min(start + block_size - 1L, n)
        block <- x[start:end, , drop = FALSE]
        centered <- sweep(block, 2, mean_vec, "-", check.margin = FALSE)
        col_ss <- col_ss + colSums(centered^2)
    }
    denom_sd <- max(1, n - 1L)
    column_sd <- sqrt(pmax(col_ss / denom_sd, 0))
    scale_vec <- if (scale) {
        out <- column_sd
        out[out == 0] <- 1
        out
    } else {
        NULL
    }

    center_vec <- if (center) mean_vec else NULL
    total_variance <- if (scale) {
        if (center) {
            as.numeric(p)
        } else {
            scale_sq <- if (!is.null(scale_vec)) scale_vec^2 else rep(1, length(col_sum_sq))
            sum(col_sum_sq / scale_sq) / n
        }
    } else {
        if (center) {
            sum(column_sd^2)
        } else {
            sum(col_sum_sq) / n
        }
    }

    random_basis <- matrix(rnorm(p * ncomp), nrow = p, ncol = ncomp)
    svd_init <- svd(random_basis, nu = ncomp, nv = 0L)
    Q <- svd_init$u[, seq_len(ncomp), drop = FALSE]

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
            end <- min(start + block_size - 1L, n)
            block <- transform_block(x[start:end, , drop = FALSE])
            if (nrow(block) == 0L) {
                next
            }
            block_times_q <- block %*% Q_mat
            accum <- accum + crossprod(block, block_times_q)
        }
        accum
    }

    convergence_flag <- FALSE
    final_delta <- NA_real_
    iterations <- 0L
    for (iter in seq_len(max_iter)) {
        iterations <- iter
        operator_Q <- apply_operator(Q)
        svd_res <- svd(operator_Q, nu = ncomp, nv = 0L)
        Q_new <- svd_res$u[, seq_len(ncomp), drop = FALSE]
        projector_diff <- Q_new %*% t(Q_new) - Q %*% t(Q)
        final_delta <- base::norm(projector_diff, type = "F")
        if (verbose) {
            message(sprintf("Iteration %d, subspace delta = %.6g", iter, final_delta))
        }
        Q <- Q_new
        if (final_delta <= tol) {
            convergence_flag <- TRUE
            break
        }
    }
    if (verbose && !convergence_flag) {
        message("Maximum iterations reached without meeting tolerance")
    }

    operator_Q <- apply_operator(Q)
    denom <- if (center) n - 1 else n
    if (denom <= 0) {
        stop("Unable to compute covariance with fewer than two observations", call. = FALSE)
    }
    small_matrix <- crossprod(Q, operator_Q) / denom
    sym_small <- (small_matrix + t(small_matrix)) / 2
    eig <- eigen(sym_small, symmetric = TRUE)
    rotation <- Q %*% eig$vectors[, seq_len(ncomp), drop = FALSE]
    sdev <- sqrt(pmax(eig$values[seq_len(ncomp)], 0))
    eigenvalues <- sdev^2

    explained <- if (total_variance > 0) eigenvalues / total_variance else rep(0, length(eigenvalues))
    cumulative <- cumsum(explained)

    colnames(rotation) <- paste0("PC", seq_len(ncomp))
    if (!is.null(colnames(x))) {
        rownames(rotation) <- colnames(x)
    }
    names(sdev) <- colnames(rotation)
    names(eigenvalues) <- colnames(rotation)
    names(explained) <- colnames(rotation)
    names(cumulative) <- colnames(rotation)

    scores <- NULL
    if (isTRUE(return_scores)) {
        scores <- matrix(0, nrow = n, ncol = ncomp)
        row_index <- 1L
        for (start in block_starts) {
            end <- min(start + block_size - 1L, n)
            block <- transform_block(x[start:end, , drop = FALSE])
            if (nrow(block) == 0L) {
                next
            }
            block_scores <- block %*% rotation
            rows <- seq.int(row_index, length.out = nrow(block))
            scores[rows, ] <- block_scores
            row_index <- row_index + nrow(block)
        }
        colnames(scores) <- colnames(rotation)
        if (!is.null(rownames(x))) {
            rownames(scores) <- rownames(x)
        }
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
    result <- new_bigpca_result(result, "spca_r")
    attr(result, "iterations") <- iterations
    attr(result, "tolerance") <- tol
    attr(result, "converged") <- isTRUE(convergence_flag)
    attr(result, "delta") <- final_delta
    result
}

#' Principal component analysis for `bigmemory::big.matrix` inputs
#' 
#' @description
#' Implements the scalable PCA (sPCA) procedure of Elgamal et al. (2015), which
#' uses block power iterations to approximate the leading principal components
#' while streaming the data in manageable chunks. The algorithm only requires
#' matrix-vector products, allowing large matrices to be processed without
#' materialising the full cross-product in memory.
#'
#' @param x A numeric matrix, data frame, [`bigmemory::big.matrix`], or an
#'   external pointer referencing a big.matrix. The input is processed in
#'   row-wise blocks so that large matrices can be analysed without creating
#'   dense copies in R memory.
#' @param ncomp Number of principal components to retain. Use `NULL` or a
#'   non-positive value to keep `min(nrow(x), ncol(x))` components.
#' @param center Logical; should column means be subtracted before performing
#'   PCA?
#' @param scale Logical; when `TRUE`, columns are scaled to unit variance after
#'   centring. Scaling requires `center = TRUE`.
#' @param block_size Number of rows to stream per block when computing column
#'   statistics and matrix-vector products.
#' @param max_iter Maximum number of block power iterations.
#' @param tol Convergence tolerance applied to the Frobenius norm of the
#'   difference between successive subspace projectors.
#' @param seed Optional integer seed used to initialise the random starting
#'   basis.
#' @param return_scores Logical; when `TRUE`, principal component scores are
#'   computed in a final streaming pass over the data.
#' @param verbose Logical; when `TRUE`, diagnostic messages describing the
#'   iteration progress are emitted.
#'
#' @return A [`bigpca`] object containing the approximate PCA solution with the
#'   same structure as [pca_bigmatrix()]. The result includes component standard
#'   deviations, rotation/loadings, optional scores, column statistics, and
#'   variance summaries. Additional metadata is stored in
#'   `attr(result, "iterations")` (number of iterations performed),
#'   `attr(result, "tolerance")` (requested tolerance), and
#'   `attr(result, "converged")` (logical convergence flag).
#'
#' @references
#' Tarek Elgamal, Maysam Yabandeh, Ashraf Aboulnaga, Waleed Mustafa, and Mohamed
#' Hefeeda (2015). *sPCA: Scalable Principal Component Analysis for Big Data on
#' Distributed Platforms*. Proceedings of the 2015 ACM SIGMOD International
#' Conference on Management of Data. <doi:10.1145/2723372.2751520>.
#'
#' @name pca_spca
#' @export
pca_spca <- function(x,
                     ncomp = NULL,
                     center = TRUE,
                     scale = FALSE,
                     block_size = 2048L,
                     max_iter = 50L,
                     tol = 1e-4,
                     seed = NULL,
                     return_scores = FALSE,
                     verbose = FALSE) {
    if (scale && !center) {
        stop("Scaling requires centring the variables", call. = FALSE)
    }

    run_spca <- function() {
        if (inherits(x, "big.matrix") || typeof(x) == "externalptr") {
            ptr <- resolve_big_pointer(x, "x")
            block_size_int <- as.integer(block_size)
            if (is.na(block_size_int) || block_size_int <= 0L) {
                stop("`block_size` must be a positive integer", call. = FALSE)
            }
            max_iter_int <- as.integer(max_iter)
            if (is.na(max_iter_int) || max_iter_int <= 0L) {
                stop("`max_iter` must be a positive integer", call. = FALSE)
            }
            ncomp_int <- if (is.null(ncomp) || ncomp <= 0L) {
                -1L
            } else {
                as.integer(ncomp)
            }
            if (is.na(ncomp_int)) {
                stop("`ncomp` must be numeric", call. = FALSE)
            }
            result <- .pca_spca_bigmatrix(
                ptr,
                center,
                scale,
                ncomp_int,
                block_size_int,
                max_iter_int,
                tol,
                return_scores,
                verbose
            )
            result <- new_bigpca_result(result, "spca_bigmemory")
            rotation <- result$rotation
            comps <- seq_len(ncol(rotation))
            colnames(rotation) <- paste0("PC", comps)
            if (inherits(x, "big.matrix")) {
                cn <- tryCatch(colnames(x), error = function(e) NULL)
                rn <- tryCatch(rownames(x), error = function(e) NULL)
                if (!is.null(cn) && length(cn) == nrow(rotation)) {
                    rownames(rotation) <- cn
                }
                if (isTRUE(return_scores) && !is.null(result$scores) &&
                    !is.null(rn) && length(rn) == nrow(result$scores)) {
                    rownames(result$scores) <- rn
                }
            }
            result$rotation <- rotation
            if (!is.null(result$sdev)) {
                names(result$sdev) <- colnames(rotation)
            }
            if (!is.null(result$eigenvalues)) {
                names(result$eigenvalues) <- colnames(rotation)
            }
            if (!is.null(result$explained_variance)) {
                names(result$explained_variance) <- colnames(rotation)
            }
            if (!is.null(result$cumulative_variance)) {
                names(result$cumulative_variance) <- colnames(rotation)
            }
            if (isTRUE(return_scores) && !is.null(result$scores)) {
                colnames(result$scores) <- colnames(rotation)
            }
            return(result)
        }
        pca_spca_R(
            x = x,
            ncomp = ncomp,
            center = center,
            scale = scale,
            block_size = block_size,
            max_iter = max_iter,
            tol = tol,
            seed = NULL,
            return_scores = return_scores,
            verbose = verbose
        )
    }

    if (is.null(seed)) {
        run_spca()
    } else {
        old_seed <- if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)) {
            get(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
        } else {
            NULL
        }
        on.exit({
            if (is.null(old_seed)) {
                if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)) {
                    rm(".Random.seed", envir = .GlobalEnv)
                }
            } else {
                assign(".Random.seed", old_seed, envir = .GlobalEnv)
            }
        }, add = TRUE)
        set.seed(seed)
        run_spca()
    }
}

#' @rdname pca_spca
#' @export
pca_spca_R <- function(x,
                       ncomp = NULL,
                       center = TRUE,
                       scale = FALSE,
                       block_size = 2048L,
                       max_iter = 50L,
                       tol = 1e-4,
                       seed = NULL,
                       return_scores = FALSE,
                       verbose = FALSE) {
    if (!is.null(seed)) {
        old_seed <- if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)) {
            get(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
        } else {
            NULL
        }
        on.exit({
            if (is.null(old_seed)) {
                if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)) {
                    rm(".Random.seed", envir = .GlobalEnv)
                }
            } else {
                assign(".Random.seed", old_seed, envir = .GlobalEnv)
            }
        }, add = TRUE)
        set.seed(seed)
    }

    if (scale && !center) {
        stop("Scaling requires centring the variables", call. = FALSE)
    }

    if (!is.matrix(x)) {
        x <- as.matrix(x)
    }
    if (!is.numeric(x)) {
        stop("`x` must contain numeric values for scalable PCA", call. = FALSE)
    }

    n <- nrow(x)
    p <- ncol(x)
    if (n == 0 || p == 0) {
        stop("`x` must contain at least one observation and one variable", call. = FALSE)
    }
    if (n < 2) {
        stop("PCA requires at least two observations", call. = FALSE)
    }

    block_size <- as.integer(block_size)
    if (is.na(block_size) || block_size <= 0L) {
        stop("`block_size` must be a positive integer", call. = FALSE)
    }
    block_size <- min(block_size, n)

    max_iter <- as.integer(max_iter)
    if (is.na(max_iter) || max_iter <= 0L) {
        stop("`max_iter` must be a positive integer", call. = FALSE)
    }

    if (is.null(ncomp) || ncomp <= 0L) {
        ncomp <- min(n, p)
    }
    ncomp <- as.integer(min(ncomp, n, p))
    if (ncomp < 1L) {
        stop("`ncomp` must select at least one component", call. = FALSE)
    }

    block_starts <- seq.int(1L, n, by = block_size)
    col_sum <- numeric(p)
    col_sum_sq <- numeric(p)
    observed <- 0L
    for (start in block_starts) {
        end <- min(start + block_size - 1L, n)
        block <- x[start:end, , drop = FALSE]
        col_sum <- col_sum + colSums(block)
        col_sum_sq <- col_sum_sq + colSums(block^2)
        observed <- observed + nrow(block)
    }
    if (observed != n) {
        n <- observed
        block_size <- min(block_size, n)
        block_starts <- seq.int(1L, n, by = block_size)
    }

    mean_vec <- col_sum / n
    col_ss <- numeric(p)
    for (start in block_starts) {
        end <- min(start + block_size - 1L, n)
        block <- x[start:end, , drop = FALSE]
        centered <- sweep(block, 2, mean_vec, "-", check.margin = FALSE)
        col_ss <- col_ss + colSums(centered^2)
    }
    denom_sd <- max(1, n - 1L)
    column_sd <- sqrt(pmax(col_ss / denom_sd, 0))
    scale_vec <- if (scale) {
        out <- column_sd
        out[out == 0] <- 1
        out
    } else {
        NULL
    }

    center_vec <- if (center) mean_vec else NULL
    total_variance <- if (scale) {
        if (center) {
            as.numeric(p)
        } else {
            scale_sq <- if (!is.null(scale_vec)) scale_vec^2 else rep(1, length(col_sum_sq))
            sum(col_sum_sq / scale_sq) / n
        }
    } else {
        if (center) {
            sum(column_sd^2)
        } else {
            sum(col_sum_sq) / n
        }
    }

    random_basis <- matrix(rnorm(p * ncomp), nrow = p, ncol = ncomp)
    svd_init <- svd(random_basis, nu = ncomp, nv = 0L)
    Q <- svd_init$u[, seq_len(ncomp), drop = FALSE]

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
            end <- min(start + block_size - 1L, n)
            block <- transform_block(x[start:end, , drop = FALSE])
            if (nrow(block) == 0L) {
                next
            }
            block_times_q <- block %*% Q_mat
            accum <- accum + crossprod(block, block_times_q)
        }
        accum
    }

    convergence_flag <- FALSE
    final_delta <- NA_real_
    iterations <- 0L
    for (iter in seq_len(max_iter)) {
        iterations <- iter
        operator_Q <- apply_operator(Q)
        svd_res <- svd(operator_Q, nu = ncomp, nv = 0L)
        Q_new <- svd_res$u[, seq_len(ncomp), drop = FALSE]
        projector_diff <- Q_new %*% t(Q_new) - Q %*% t(Q)
        final_delta <- base::norm(projector_diff, type = "F")
        if (verbose) {
            message(sprintf("Iteration %d, subspace delta = %.6g", iter, final_delta))
        }
        Q <- Q_new
        if (final_delta <= tol) {
            convergence_flag <- TRUE
            break
        }
    }
    if (verbose && !convergence_flag) {
        message("Maximum iterations reached without meeting tolerance")
    }

    operator_Q <- apply_operator(Q)
    denom <- if (center) n - 1 else n
    if (denom <= 0) {
        stop("Unable to compute covariance with fewer than two observations", call. = FALSE)
    }
    small_matrix <- crossprod(Q, operator_Q) / denom
    sym_small <- (small_matrix + t(small_matrix)) / 2
    eig <- eigen(sym_small, symmetric = TRUE)
    rotation <- Q %*% eig$vectors[, seq_len(ncomp), drop = FALSE]
    sdev <- sqrt(pmax(eig$values[seq_len(ncomp)], 0))
    eigenvalues <- sdev^2

    explained <- if (total_variance > 0) eigenvalues / total_variance else rep(0, length(eigenvalues))
    cumulative <- cumsum(explained)

    colnames(rotation) <- paste0("PC", seq_len(ncomp))
    if (!is.null(colnames(x))) {
        rownames(rotation) <- colnames(x)
    }
    names(sdev) <- colnames(rotation)
    names(eigenvalues) <- colnames(rotation)
    names(explained) <- colnames(rotation)
    names(cumulative) <- colnames(rotation)

    scores <- NULL
    if (isTRUE(return_scores)) {
        scores <- matrix(0, nrow = n, ncol = ncomp)
        row_index <- 1L
        for (start in block_starts) {
            end <- min(start + block_size - 1L, n)
            block <- transform_block(x[start:end, , drop = FALSE])
            if (nrow(block) == 0L) {
                next
            }
            block_scores <- block %*% rotation
            rows <- seq.int(row_index, length.out = nrow(block))
            scores[rows, ] <- block_scores
            row_index <- row_index + nrow(block)
        }
        colnames(scores) <- colnames(rotation)
        if (!is.null(rownames(x))) {
            rownames(scores) <- rownames(x)
        }
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
    result <- new_bigpca_result(result, "spca_r")
    attr(result, "iterations") <- iterations
    attr(result, "tolerance") <- tol
    attr(result, "converged") <- isTRUE(convergence_flag)
    attr(result, "delta") <- final_delta
    result
}

#' Principal component analysis for `bigmemory::big.matrix` inputs
#' 
#' @description
#' Perform principal component analysis (PCA) directly on a
#' [`bigmemory::big.matrix`] without copying the data into R memory. The
#' exported helpers mirror the structure of base R's `prcomp()` while avoiding
#' the need to materialise large matrices.
#' 
#' @name pca_bigmatrix
#' @param xpMat Either a [`bigmemory::big.matrix`] or an external pointer such
#'   as `mat@address` that references the source `big.matrix`.
#' @param center Logical; should column means be subtracted before PCA?
#' @param scale Logical; should variables be scaled to unit variance? Scaling
#'   requires `center = TRUE`.
#' @param ncomp Number of components to retain. Use a non-positive value to keep
#'   all components returned by the decomposition.
#' @param block_size Number of rows to process per block when streaming data
#'   through BLAS kernels. Larger values improve throughput at the cost of
#'   additional memory.
#'
#' @return For `pca_bigmatrix()`, a `bigpca` object mirroring a `prcomp` result
#'   with elements `sdev`, `rotation`, optional `center` and `scale` vectors,
#'   `column_sd`, `eigenvalues`, `explained_variance`, `cumulative_variance`, and
#'   the sample covariance matrix. The object participates in S3 generics such as
#'   [summary()] and [plot()].
#'
#' @seealso [bigpca], [pca_scores_bigmatrix()], [pca_variable_loadings()],
#'   [pca_variable_correlations()], [pca_variable_contributions()], and the
#'   streaming variants [pca_stream_bigmatrix()] and companions.
#' @examplesIf requireNamespace("bigmemory", quietly = TRUE)
#' set.seed(123)
#' mat <- bigmemory::as.big.matrix(matrix(rnorm(40), nrow = 10))
#' pca <- pca_bigmatrix(mat, center = TRUE, scale = TRUE, ncomp = 3)
#' scores <- pca_scores_bigmatrix(mat, pca$rotation, pca$center, pca$scale, ncomp = 3)
#' loadings <- pca_variable_loadings(pca$rotation, pca$sdev)
#' correlations <- pca_variable_correlations(pca$rotation, pca$sdev, pca$column_sd)
#' contributions <- pca_variable_contributions(loadings)
#' list(scores = scores, loadings = loadings, correlations = correlations,
#'      contributions = contributions)
#' @param arg Character string naming the argument being validated. Used to
#'   construct informative error messages.
#' @param allow_null Logical flag indicating whether `NULL` is accepted for the
#'   argument. When `TRUE`, a `NULL` input is returned unchanged.
resolve_big_pointer <- function(x, arg, allow_null = FALSE) {
    if (allow_null && is.null(x)) {
        return(NULL)
    }
    if (methods::is(x, "big.matrix")) {
        addr <- methods::slot(x, "address")
        if (!identical(typeof(addr), "externalptr")) {
            stop(sprintf("`%s@address` must be an external pointer", arg), call. = FALSE)
        }
        return(addr)
    }
    if (typeof(x) == "externalptr") {
        return(x)
    }
    stop(sprintf("`%s` must be a bigmemory::big.matrix or external pointer", arg), call. = FALSE)
}

#' Singular value decomposition for `bigmemory::big.matrix` inputs
#'
#' @description
#' Compute the singular value decomposition (SVD) of a
#' [`bigmemory::big.matrix`] without materialising it as a base R matrix.
#' Blocks of rows are streamed through BLAS before LAPACK is invoked so that
#' even moderately large matrices can be decomposed efficiently.
#'
#' @param xpMat Either a [`bigmemory::big.matrix`] or an external pointer such
#'   as `mat@address` that references the source `big.matrix`.
#' @param nu Number of left singular vectors to return. Use a negative value to
#'   request the default of `min(nrow, ncol)` vectors and zero to skip
#'   returning `u` entirely.
#' @param nv Number of right singular vectors to return. Use a negative value to
#'   request the default of `min(nrow, ncol)` vectors and zero to skip
#'   returning `v` entirely.
#' @param block_size Number of rows to process per block when streaming data
#'   into BLAS kernels. Larger values can improve throughput at the cost of
#'   additional temporary memory.
#' @param method LAPACK backend used to compute the decomposition. The default
#'   uses the divide-and-conquer routine `dgesdd` and falls back to `dgesvd`
#'   when required.
#'
#' @return A list with components `u`, `d`, and `v` analogous to base R's
#'   [svd()] output. When `nu` or `nv` are zero the corresponding matrix has
#'   zero columns.
#'
#' @examplesIf requireNamespace("bigmemory", quietly = TRUE)
#' set.seed(42)
#' mat <- bigmemory::as.big.matrix(matrix(rnorm(20), nrow = 5))
#' svd_res <- svd_bigmatrix(mat, nu = 2, nv = 2)
#' svd_res$d
#'
#' @export
svd_bigmatrix <- function(xpMat, nu = -1L, nv = -1L, block_size = 1024L, method = c("dgesdd", "dgesvd")) {
    ptr <- resolve_big_pointer(xpMat, "xpMat")
    method <- match.arg(method)
    .svd_bigmatrix(ptr, as.integer(nu), as.integer(nv), as.integer(block_size), method)
}
#' @export
pca_bigmatrix <- function(xpMat, center = TRUE, scale = FALSE, ncomp = -1L, block_size = 1024L) {
    ptr <- resolve_big_pointer(xpMat, "xpMat")
    result <- .pca_bigmatrix(ptr, center, scale, ncomp, block_size)
    new_bigpca_result(result, "bigmemory")
}

#' @describeIn pca_bigmatrix Project observations into principal component
#'   space while streaming from a `big.matrix`.
# #' @inheritParams pca_bigmatrix
#' @param rotation A rotation matrix such as the `rotation` element returned by
#'   [pca_bigmatrix()].
#' @param center For `pca_scores_bigmatrix()`, a numeric vector of column means
#'   (optional).
#' @param scale For `pca_scores_bigmatrix()`, a numeric vector of scaling
#'   factors (optional). Scaling requires `center` to be supplied.
#' @return A numeric matrix of scores with rows corresponding to observations
#'   and columns to retained components.
#' @export
pca_scores_bigmatrix <- function(xpMat, rotation, center, scale, ncomp = -1L, block_size = 1024L) {
    ptr <- resolve_big_pointer(xpMat, "xpMat")
    if (is.null(center)) {
        center <- numeric()
    }
    if (is.null(scale)) {
        scale <- numeric()
    }
    .pca_scores_bigmatrix(ptr, rotation, center, scale, ncomp, block_size)
}

#' @describeIn pca_bigmatrix Compute variable loadings (covariances between
#'   original variables and components).
#' @param sdev A numeric vector of component standard deviations, typically the
#'   `sdev` element from [pca_bigmatrix()].
#' @return A numeric matrix containing variable loadings for each component.
#' @export
pca_variable_loadings <- function(rotation, sdev) {
    .pca_variable_loadings(rotation, sdev)
}

#' @describeIn pca_bigmatrix Compute variable-component correlations given
#'   column standard deviations.
#' @param column_sd A numeric vector with the marginal standard deviation of
#'   each original variable.
#' @return A numeric matrix of correlations between variables and components.
#' @export
pca_variable_correlations <- function(rotation, sdev, column_sd) {
    .pca_variable_correlations(rotation, sdev, column_sd)
}

#' @describeIn pca_bigmatrix Derive the relative contribution of each variable
#'   to the retained components.
#' @param loadings A numeric matrix such as the result of
#'   [pca_variable_loadings()].
#' @param scores For `pca_individual_contributions()` and
#'   `pca_individual_cos2()`, a numeric matrix of component scores where rows
#'   correspond to observations and columns to components.
#' @param total_weight Optional positive scalar giving the effective number of
#'   observations to use when computing contributions. Defaults to the number of
#'   rows in `scores`.
#' @return A numeric matrix where each entry represents the contribution of a
#'   variable to a component.
#' @export
pca_variable_contributions <- function(loadings) {
    .pca_variable_contributions(loadings)
}

#' @describeIn pca_bigmatrix Compute the relative contribution of individual
#'   observations to each component.
#' @export
pca_individual_contributions <- function(scores, sdev, total_weight = NA_real_) {
    .pca_individual_contributions(scores, sdev, total_weight)
}

#' @describeIn pca_bigmatrix Compute squared cosine values measuring the quality
#'   of representation for individual observations.
#' @export
pca_individual_cos2 <- function(scores) {
    .pca_individual_cos2(scores)
}

#' @describeIn pca_bigmatrix Compute squared cosine values measuring the quality
#'   of representation for variables.
#' @param correlations For `pca_variable_cos2()`, a numeric matrix of
#'   correlations between variables and components.
#' @export
pca_variable_cos2 <- function(correlations) {
    .pca_variable_cos2(correlations)
}

#' Supplementary individual diagnostics
#'
#' @description Compute principal component scores and quality metrics for
#'   supplementary individuals (rows) projected into an existing PCA solution.
#'
#' @param data Matrix-like object whose rows correspond to supplementary
#'   individuals and columns to the original variables.
#' @param rotation Rotation matrix from the PCA model (e.g. the `rotation`
#'   element of a [`bigpca`] result).
#' @param sdev Numeric vector of component standard deviations associated with
#'   `rotation`.
#' @param center Optional numeric vector giving the centring applied to each
#'   variable when fitting the PCA. Defaults to zero centring.
#' @param scale Optional numeric vector describing the scaling applied to each
#'   variable when fitting the PCA. When `NULL`, no scaling is applied.
#' @param total_weight Optional positive scalar passed to
#'   [pca_individual_contributions()] when computing contributions. When left as
#'   `NA` (the default), the resulting contributions for each component are
#'   normalised to sum to one across supplementary individuals. Supplying a
#'   value bypasses this normalisation and delegates the scaling to
#'   `pca_individual_contributions()`.
#' @return A list with elements `scores`, `contributions`, and `cos2`.
#' @export
pca_supplementary_individuals <- function(data, rotation, sdev, center = NULL, scale = NULL,
                                          total_weight = NA_real_) {
    mat <- as.matrix(data)
    rotation <- as.matrix(rotation)
    if (!is.numeric(mat)) {
        stop("`data` must be coercible to a numeric matrix", call. = FALSE)
    }
    if (!is.numeric(rotation)) {
        stop("`rotation` must be a numeric matrix", call. = FALSE)
    }
    if (ncol(mat) != nrow(rotation)) {
        stop("Number of columns in `data` must match rows of `rotation`", call. = FALSE)
    }
    k <- ncol(rotation)
    if (length(sdev) < k) {
        stop("`sdev` must contain at least as many elements as there are components", call. = FALSE)
    }

    center_vec <- resolve_center_vector(center, ncol(mat))
    scale_vec <- resolve_scale_vector(scale, ncol(mat))

    centered <- sweep(mat, 2, center_vec, "-", check.margin = FALSE)
    scaled <- sweep(centered, 2, scale_vec, "/", check.margin = FALSE)

    scores <- scaled %*% rotation
    if (!is.null(rownames(mat))) {
        rownames(scores) <- rownames(mat)
    }
    if (!is.null(colnames(rotation))) {
        colnames(scores) <- colnames(rotation)
    }

    contributions <- if (isTRUE(is.na(total_weight))) {
        contrib <- pca_individual_contributions(scores, sdev)
        if (is.null(dim(contrib))) {
            contrib <- matrix(contrib, nrow = nrow(scores), ncol = ncol(scores))
        }
        totals <- colSums(contrib)
        if (length(totals)) {
            for (comp in seq_along(totals)) {
                total <- totals[comp]
                if (total > 0) {
                    contrib[, comp] <- contrib[, comp] / total
                }
            }
        }
        contrib
    } else {
        contrib <- pca_individual_contributions(scores, sdev, total_weight = total_weight)
        if (is.null(dim(contrib))) {
            contrib <- matrix(contrib, nrow = nrow(scores), ncol = ncol(scores))
        }
        totals <- colSums(contrib)
        if (length(totals)) {
            for (comp in seq_along(totals)) {
                total <- totals[comp]
                if (total > 0) {
                    contrib[, comp] <- contrib[, comp] * (total_weight / total)
                }
            }
        }
        contrib
    }
    cos2 <- pca_individual_cos2(scores)

    if (!is.null(rownames(scores))) {
        rownames(contributions) <- rownames(scores)
        rownames(cos2) <- rownames(scores)
    }
    if (!is.null(colnames(scores))) {
        colnames(contributions) <- colnames(scores)
        colnames(cos2) <- colnames(scores)
    }

    list(scores = scores, contributions = contributions, cos2 = cos2)
}

#' Supplementary variable diagnostics
#'
#' @description Compute loadings, correlations, contributions, and cos^2 values
#'   for supplementary variables (columns) given component scores for the active
#'   individuals.
#'
#' @param data Matrix-like object whose columns correspond to supplementary
#'   variables measured on the active individuals.
#' @param scores Numeric matrix of component scores for the active individuals.
#' @param sdev Numeric vector of component standard deviations associated with
#'   `scores`.
#' @param center Optional numeric vector specifying the centring to apply to each
#'   supplementary variable. When `NULL`, column means of `data` are used.
#' @return A list with elements `loadings`, `correlations`, `contributions`, and
#'   `cos2`.
#' @export
pca_supplementary_variables <- function(data, scores, sdev, center = NULL) {
    mat <- as.matrix(data)
    score_mat <- as.matrix(scores)
    if (!is.numeric(mat)) {
        stop("`data` must be coercible to a numeric matrix", call. = FALSE)
    }
    if (!is.numeric(score_mat)) {
        stop("`scores` must be a numeric matrix", call. = FALSE)
    }
    if (nrow(mat) != nrow(score_mat)) {
        stop("Supplementary variables must be observed on the same individuals as `scores`", call. = FALSE)
    }
    k <- ncol(score_mat)
    if (length(sdev) < k) {
        stop("`sdev` must contain at least as many elements as there are components", call. = FALSE)
    }

    center_vec <- resolve_center_vector(center, ncol(mat), default = NA_real_)
    if (anyNA(center_vec)) {
        center_vec <- colMeans(mat)
    }

    centered <- sweep(mat, 2, center_vec, "-", check.margin = FALSE)
    score_center <- colMeans(score_mat)
    scores_centered <- sweep(score_mat, 2, score_center, "-", check.margin = FALSE)

    cross <- crossprod(centered, scores_centered)
    var_norm <- sqrt(colSums(centered^2))
    score_norm <- sqrt(colSums(scores_centered^2))
    denom <- outer(var_norm, score_norm)
    correlations <- cross
    valid <- denom > .Machine$double.eps
    correlations[valid] <- correlations[valid] / denom[valid]
    correlations[!valid] <- 0

    loadings <- correlations
    loadings <- sweep(loadings, 2, sdev[seq_len(k)], `*`)

    if (!is.null(colnames(score_mat))) {
        colnames(loadings) <- colnames(score_mat)
        colnames(correlations) <- colnames(score_mat)
    }
    if (!is.null(colnames(mat))) {
        rownames(loadings) <- colnames(mat)
        rownames(correlations) <- colnames(mat)
    }

    contributions <- pca_variable_contributions(loadings)
    cos2 <- pca_variable_cos2(correlations)

    if (!is.null(rownames(loadings))) {
        rownames(contributions) <- rownames(loadings)
        rownames(cos2) <- rownames(loadings)
    }
    if (!is.null(colnames(loadings))) {
        colnames(contributions) <- colnames(loadings)
        colnames(cos2) <- colnames(loadings)
    }

    list(
        loadings = loadings,
        correlations = correlations,
        contributions = contributions,
        cos2 = cos2
    )
}

resolve_center_vector <- function(center, n, default = 0) {
    if (is.null(center)) {
        return(rep(default, n))
    }
    center_vec <- as.numeric(center)
    if (!length(center_vec)) {
        return(rep(default, n))
    }
    if (length(center_vec) == 1L) {
        center_vec <- rep(center_vec, n)
    }
    if (length(center_vec) != n) {
        stop("Length of `center` must match the number of variables", call. = FALSE)
    }
    center_vec
}

resolve_scale_vector <- function(scale, n) {
    if (is.null(scale)) {
        return(rep(1, n))
    }
    scale_vec <- as.numeric(scale)
    if (!length(scale_vec)) {
        return(rep(1, n))
    }
    if (length(scale_vec) == 1L) {
        scale_vec <- rep(scale_vec, n)
    }
    if (length(scale_vec) != n) {
        stop("Length of `scale` must match the number of variables", call. = FALSE)
    }
    if (any(scale_vec == 0)) {
        stop("`scale` entries must be non-zero", call. = FALSE)
    }
    scale_vec
}

#' Streaming big.matrix PCA helpers
#'
#' @description
#' Variants of the PCA helpers that stream results directly into
#' `bigmemory::big.matrix` objects, enabling file-backed workflows without
#' materialising dense R matrices.
#'
#' @inheritParams pca_bigmatrix
#' @inheritParams pca_spca
#' @param xpRotation Optionally, either a `bigmemory::big.matrix` or external
#'   pointer referencing a destination big.matrix that receives the rotation
#'   matrix.
#' @return For `pca_stream_bigmatrix()`, the same `bigpca` object as
#'   [pca_bigmatrix()] with the
#'   addition of a `rotation_stream_bigmatrix` element referencing the populated
#'   `big.matrix` when `xpRotation` is supplied. For
#'   `pca_spca_stream_bigmatrix()`, the same scalable PCA structure as
#'   [pca_spca()] with the optional pointer populated when provided.
#' @name pca_stream_bigmatrix
#' @export
pca_spca_stream_bigmatrix <- function(
        xpMat,
        xpRotation = NULL,
        center = TRUE,
        scale = FALSE,
        ncomp = -1L,
        block_size = 2048L,
        max_iter = 50L,
        tol = 1e-4,
        seed = NULL,
        return_scores = FALSE,
        verbose = FALSE) {
    mat_ptr <- resolve_big_pointer(xpMat, "xpMat")
    rot_ptr <- resolve_big_pointer(xpRotation, "xpRotation", allow_null = TRUE)
    block_size <- as.integer(block_size)
    if (is.na(block_size) || block_size <= 0L) {
        stop("`block_size` must be a positive integer", call. = FALSE)
    }
    max_iter <- as.integer(max_iter)
    if (is.na(max_iter) || max_iter <= 0L) {
        stop("`max_iter` must be a positive integer", call. = FALSE)
    }
    if (!is.null(seed)) {
        old_seed <- if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)) {
            get(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
        } else {
            NULL
        }
        on.exit({
            if (is.null(old_seed)) {
                if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)) {
                    rm(".Random.seed", envir = .GlobalEnv)
                }
            } else {
                assign(".Random.seed", old_seed, envir = .GlobalEnv)
            }
        }, add = TRUE)
        set.seed(seed)
    }
    result <- .pca_spca_stream_bigmatrix(mat_ptr, rot_ptr, center, scale, ncomp, block_size, max_iter, tol, return_scores, verbose)
    result <- new_bigpca_result(result, "spca_stream_bigmatrix")
    rotation <- result$rotation
    comps <- seq_len(ncol(rotation))
    colnames(rotation) <- paste0("PC", comps)
    if (methods::is(xpMat, "big.matrix")) {
        cn <- tryCatch(colnames(xpMat), error = function(e) NULL)
        rn <- tryCatch(rownames(xpMat), error = function(e) NULL)
        if (!is.null(cn) && length(cn) == nrow(rotation)) {
            rownames(rotation) <- cn
        }
        if (isTRUE(return_scores) && !is.null(result$scores) && !is.null(rn) && length(rn) == nrow(result$scores)) {
            rownames(result$scores) <- rn
        }
    }
    result$rotation <- rotation
    if (!is.null(result$sdev)) {
        names(result$sdev) <- colnames(rotation)
    }
    if (!is.null(result$eigenvalues)) {
        names(result$eigenvalues) <- colnames(rotation)
    }
    if (!is.null(result$explained_variance)) {
        names(result$explained_variance) <- colnames(rotation)
    }
    if (!is.null(result$cumulative_variance)) {
        names(result$cumulative_variance) <- colnames(rotation)
    }
    if (isTRUE(return_scores) && !is.null(result$scores)) {
        colnames(result$scores) <- colnames(rotation)
    }
    result
}

#' @export
pca_stream_bigmatrix <- function(
        xpMat,
        xpRotation = NULL,
        center = TRUE,
        scale = FALSE,
        ncomp = -1L,
        block_size = 1024L) {
    mat_ptr <- resolve_big_pointer(xpMat, "xpMat")
    rot_ptr <- resolve_big_pointer(xpRotation, "xpRotation", allow_null = TRUE)
    result <- .pca_stream_bigmatrix(mat_ptr, rot_ptr, center, scale, ncomp, block_size)
    new_bigpca_result(result, "stream_bigmatrix")
}

#' @describeIn pca_stream_bigmatrix Stream PCA scores into a destination
#'   big.matrix.
#' @param xpDest A `bigmemory::big.matrix` or external pointer referencing a
#'   writable destination `big.matrix` that will receive the projected scores.
#' @export
pca_scores_stream_bigmatrix <- function(
        xpMat,
        xpDest,
        rotation,
        center,
        scale,
        ncomp = -1L,
        block_size = 1024L) {
    mat_ptr <- resolve_big_pointer(xpMat, "xpMat")
    dest_ptr <- resolve_big_pointer(xpDest, "xpDest")
    if (is.null(center)) {
        center <- numeric()
    }
    if (is.null(scale)) {
        scale <- numeric()
    }
    .pca_scores_stream_bigmatrix(mat_ptr, dest_ptr, rotation, center, scale, ncomp, block_size)
}

#' @describeIn pca_stream_bigmatrix Populate big.matrix objects with derived variable
#'   diagnostics.
#' @param xpLoadings For `pca_variable_contributions_stream_bigmatrix()`, the
#'   loadings matrix supplied as a `big.matrix` or external pointer.
#' @param xpDest Either a `big.matrix` or external pointer referencing the
#'   destination `big.matrix` that stores the computed quantity.
#' @return The external pointer supplied in `xpDest`, invisibly.
#' @examplesIf requireNamespace("bigmemory", quietly = TRUE)
#' set.seed(456)
#' mat <- bigmemory::as.big.matrix(matrix(rnorm(30), nrow = 6))
#' ncomp <- 2
#' rotation_store <- bigmemory::big.matrix(ncol(mat), ncomp, type = "double")
#' pca_stream <- pca_stream_bigmatrix(mat, xpRotation = rotation_store, ncomp = ncomp)
#' score_store <- bigmemory::big.matrix(nrow(mat), ncomp, type = "double")
#' pca_scores_stream_bigmatrix(
#'     mat,
#'     score_store,
#'     pca_stream$rotation,
#'     pca_stream$center,
#'     pca_stream$scale,
#'     ncomp = ncomp
#' )
#' loadings_store <- bigmemory::big.matrix(ncol(mat), ncomp, type = "double")
#' pca_variable_loadings_stream_bigmatrix(
#'     pca_stream$rotation_stream_bigmatrix,
#'     pca_stream$sdev,
#'     loadings_store
#' )
#' correlation_store <- bigmemory::big.matrix(ncol(mat), ncomp, type = "double")
#' pca_variable_correlations_stream_bigmatrix(
#'     pca_stream$rotation_stream_bigmatrix,
#'     pca_stream$sdev,
#'     pca_stream$column_sd,
#'     correlation_store
#' )
#' contribution_store <- bigmemory::big.matrix(ncol(mat), ncomp, type = "double")
#' pca_variable_contributions_stream_bigmatrix(
#'     loadings_store,
#'     contribution_store
#' )
#' @export
pca_variable_loadings_stream_bigmatrix <- function(
        xpRotation,
        sdev,
        xpDest) {
    rot_ptr <- resolve_big_pointer(xpRotation, "xpRotation")
    dest_ptr <- resolve_big_pointer(xpDest, "xpDest")
    .pca_variable_loadings_stream_bigmatrix(rot_ptr, sdev, dest_ptr)
}

#' @describeIn pca_stream_bigmatrix Stream variable correlations into a
#'   destination big.matrix.
#' @param xpRotation For `pca_variable_correlations_stream_bigmatrix()`, a
#'   `bigmemory::big.matrix` or external pointer containing the rotation matrix
#'   to stream from.
#' @param column_sd A numeric vector of variable standard deviations used to
#'   scale the correlations.
#' @export
pca_variable_correlations_stream_bigmatrix <- function(
        xpRotation,
        sdev,
        column_sd,
        xpDest) {
    rot_ptr <- resolve_big_pointer(xpRotation, "xpRotation")
    dest_ptr <- resolve_big_pointer(xpDest, "xpDest")
    .pca_variable_correlations_stream_bigmatrix(rot_ptr, sdev, column_sd, dest_ptr)
}

#' @describeIn pca_stream_bigmatrix Stream variable contributions into a
#'   destination big.matrix.
#' @export
pca_variable_contributions_stream_bigmatrix <- function(
        xpLoadings,
        xpDest) {
    load_ptr <- resolve_big_pointer(xpLoadings, "xpLoadings")
    dest_ptr <- resolve_big_pointer(xpDest, "xpDest")
    .pca_variable_contributions_stream_bigmatrix(load_ptr, dest_ptr)
}
