#' BigPCA result objects
#'
#' @description
#' Results returned by [pca_bigmatrix()], [pca_stream_bigmatrix()], and
#' [pca_robust()] inherit from the `bigpca` class. The objects store the
#' component standard deviations, rotation/loadings, and optional scores while
#' recording which computational backend produced them. Standard S3 generics
#' such as [summary()] and [plot()] are implemented for convenience.
#'   
#' `bigpca` objects are lists produced by [pca_bigmatrix()],
#' [pca_stream_bigmatrix()], [pca_robust()], and related helpers. They mirror
#' the structure of base R's `prcomp()` outputs while tracking additional
#' metadata for large-scale and streaming computations.
#' 
#' #' @seealso [pca_bigmatrix()], [pca_stream_bigmatrix()], [pca_robust()],
#'   [pca_plot_scree()], [pca_plot_scores()], [pca_plot_contributions()],
#'   [pca_plot_correlation_circle()], and [pca_plot_biplot()].
#'
#' @section Components:
#' \describe{
#'   \item{`sdev`}{Numeric vector of component standard deviations.}
#'   \item{`rotation`}{Numeric matrix whose columns contain the variable
#'   loadings (principal axes).}
#'   \item{`center`, `scale`}{Optional numeric vectors describing the centring
#'   and scaling applied to each variable when fitting the model.}
#'   \item{`scores`}{Optional numeric matrix of principal component scores when
#'   computed alongside the decomposition.}
#'   \item{`column_sd`}{Numeric vector of marginal standard deviations for each
#'   input variable.}
#'   \item{`eigenvalues`}{Numeric vector of eigenvalues associated with the
#'   retained components.}
#'   \item{`explained_variance`, `cumulative_variance`}{Numeric vectors
#'   summarising the fraction of variance explained by individual components
#'   and the corresponding cumulative totals.}
#'   \item{`covariance`}{Sample covariance matrix used to derive the
#'   components.}
#'   \item{`nobs`}{Number of observations used in the decomposition.}
#' }
#'
#' The class also records the computation backend via
#' `attr(x, "backend")`, enabling downstream methods to adjust their
#' behaviour for streamed or robust results.
#'
#' @seealso [pca_bigmatrix()], [pca_stream_bigmatrix()], [summary.bigpca()],
#'   [print.summary.bigpca()], [plot.bigpca()]
#' @family bigpca
#' @name bigpca
#' @rdname bigpca
NULL

#' Internal constructors and S3 methods for bigPCAcpp results
#'
#' These helpers provide a light-weight S3 layer around the PCA outputs so
#' users can interact with them through familiar generics such as
#' [summary()] and [plot()].
#'
#' @keywords internal
new_bigpca_result <- function(result, backend) {
    stopifnot(is.list(result))
    if (!is.character(backend) || length(backend) != 1L) {
        stop("`backend` must be a single character string")
    }
    attr(result, "backend") <- backend
    class(result) <- unique(c(sprintf("bigpca_%s", backend), "bigpca", class(result)))
    result
}

#' @describeIn pca_bigmatrix Summarise the component importance metrics for a
#'   [`bigpca`] result.
#' @param object A [`bigpca`] object created by [pca_bigmatrix()],
#'   [pca_stream_bigmatrix()], or related helpers.
#' @param ... Additional arguments passed to methods. Currently unused.
#' @return For `summary.bigpca()`, a [`summary.bigpca`] object containing
#'   component importance measures.
#' @seealso [bigpca]
#' @export
summary.bigpca <- function(object, ...) {
    if (is.null(object$sdev)) {
        stop("`object$sdev` is required to summarise a PCA result")
    }
    sdev <- as.numeric(object$sdev)
    if (!length(sdev)) {
        stop("`object$sdev` must contain at least one component")
    }
    explained <- object$explained_variance
    if (is.null(explained)) {
        total <- sum(sdev^2)
        explained <- if (total > 0) (sdev^2) / total else rep(0, length(sdev))
    } else {
        explained <- as.numeric(explained)
    }
    cumulative <- object$cumulative_variance
    if (is.null(cumulative)) {
        cumulative <- cumsum(explained)
    } else {
        cumulative <- as.numeric(cumulative)
    }
    importance <- rbind(
        `Standard deviation` = sdev,
        `Proportion of Variance` = explained,
        `Cumulative Proportion` = cumulative
    )
    backend <- attr(object, "backend", exact = TRUE)
    if (is.null(backend) && !is.null(object$backend)) {
        backend <- object$backend
    }
    result <- list(
        importance = importance,
        backend = backend %||% "unknown",
        nobs = object$nobs
    )
    class(result) <- "summary.bigpca"
    result
}

#' @describeIn pca_bigmatrix Print the component importance summary produced by
#'   [summary.bigpca()].
#' @param x A [`summary.bigpca`] object.
#' @param digits Number of significant digits to display when printing
#'   importance metrics.
#' @seealso [bigpca]
#' @export
print.summary.bigpca <- function(x, digits = max(3, getOption("digits") - 3), ...) {
    cat("Big PCA summary", if (!is.null(x$backend)) sprintf("(%s)", x$backend) else "", "\n", sep = "")
    if (!is.null(x$nobs)) {
        cat("Observations:", format(x$nobs, big.mark = ","), "\n")
    }
    importance <- x$importance
    if (!is.null(importance)) {
        print(round(importance, digits = digits), quote = FALSE)
    }
    invisible(x)
}

#' @describeIn pca_bigmatrix Visualise PCA diagnostics such as scree, correlation
#'   circle, contribution, and biplot displays.
#' @param y Currently unused.
#' @param type The plot to draw. Options include "scree" (variance explained),
#'   "contributions" (top contributing variables), "correlation_circle" (variable
#'   correlations with selected components), and "biplot" (joint display of
#'   scores and loadings).
#' @param max_components Maximum number of components to display in scree plots.
#' @param component Component index to highlight when drawing contribution plots.
#' @param top_n Number of variables to display in contribution plots.
#' @param components Length-two integer vector selecting the components for
#'   correlation circle and biplot views.
#' @param data Optional data source (matrix, data frame, `bigmemory::big.matrix`,
#'   or external pointer) used to compute scores for biplots when
#'   `x$scores` is unavailable.
#' @param draw Logical; if `FALSE`, return the data prepared for the selected
#'   plot instead of drawing it.
#' @param ... Additional arguments passed to plotting helpers.
#' @seealso [bigpca]
#' @export
plot.bigpca <- function(x, y, type = c("scree", "contributions", "correlation_circle", "biplot"),
                        max_components = 25L, component = 1L, top_n = 20L,
                        components = c(1L, 2L), data = NULL, draw = TRUE, ...) {
    type <- match.arg(type)
    if (identical(type, "scree")) {
        return(invisible(pca_plot_scree(x, max_components = max_components, draw = draw, ...)))
    }

    if (identical(type, "contributions")) {
        loadings <- pca_variable_loadings(x$rotation, x$sdev)
        contributions <- pca_variable_contributions(loadings)
        return(invisible(pca_plot_contributions(contributions, component = component,
                                                top_n = top_n, draw = draw, ...)))
    }

    if (identical(type, "correlation_circle")) {
        correlations <- pca_variable_correlations(x$rotation, x$sdev, x$column_sd, x$scale)
        return(invisible(pca_plot_correlation_circle(correlations, components = components,
                                                     draw = draw, ...)))
    }

    loadings <- pca_variable_loadings(x$rotation, x$sdev)
    scores <- extract_biplot_scores(x, data, components)
    invisible(pca_plot_biplot(scores, loadings, components = components, draw = draw, ...))
}

extract_biplot_scores <- function(result, data, components) {
    if (length(components) != 2L) {
        stop("`components` must contain exactly two component indices")
    }
    max_comp <- max(components)
    rotation <- as.matrix(result$rotation)
    if (ncol(rotation) < max_comp) {
        stop("Requested components are not available in the PCA rotation")
    }
    rotation <- rotation[, seq_len(max_comp), drop = FALSE]

    scores_source <- result$scores
    if (!is.null(scores_source)) {
        scores_mat <- as_matrix(scores_source)
        if (ncol(scores_mat) < max_comp) {
            stop("Stored scores do not contain the requested components")
        }
        return(scores_mat[, components, drop = FALSE])
    }

    if (is.null(data)) {
        stop("`data` must be supplied when `x$scores` are unavailable for biplots")
    }

    center <- if (is.null(result$center)) numeric() else as.numeric(result$center)
    scale <- if (is.null(result$scale)) numeric() else as.numeric(result$scale)

    if (methods::is(data, "big.matrix") || identical(typeof(data), "externalptr")) {
        scores_full <- pca_scores_bigmatrix(data, rotation, center, scale, ncomp = max_comp)
        scores_mat <- as_matrix(scores_full)
    } else {
        mat <- as.matrix(data)
        if (ncol(mat) != nrow(rotation)) {
            stop("`data` must have the same number of columns as the PCA rotation")
        }
        if (length(center)) {
            mat <- sweep(mat, 2, center, "-")
        }
        if (length(scale)) {
            safe_scale <- ifelse(scale == 0, 1, scale)
            mat <- sweep(mat, 2, safe_scale, "/")
        }
        scores_mat <- mat %*% rotation
    }

    scores_mat[, components, drop = FALSE]
}

as_matrix <- function(x) {
    if (methods::is(x, "big.matrix")) {
        return(x[, , drop = FALSE])
    }
    if (is.matrix(x)) {
        return(x)
    }
    if (is.data.frame(x)) {
        return(as.matrix(x))
    }
    if (inherits(x, "big.matrix.descriptor")) {
        stop("`big.matrix.descriptor` objects are not supported directly; pass the backing big.matrix instead")
    }
    if (inherits(x, "Matrix")) {
        return(as.matrix(x))
    }
    as.matrix(x)
}

`%||%` <- function(x, y) {
    if (is.null(x)) y else x
}
