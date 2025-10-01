#' Plot PCA diagnostics for big data workflows
#'
#' These helpers visualise the results returned by [pca_bigmatrix()] and its
#' companions without requiring users to materialise dense intermediate
#' structures. Each plotting function optionally samples the inputs so the
#' default output remains responsive even when the underlying big matrix spans
#' millions of observations.
#'
#' @name pca_plots
#' @seealso [pca_bigmatrix()], [pca_variable_loadings()],
#'   [pca_variable_contributions()]
NULL

#' Scree plot for principal component importance
#'
#' @description
#' Displays the proportion of variance explained by the leading principal
#' components. The function caps the number of displayed components to keep the
#' visualization legible on very high-dimensional problems.
#'
#' @param pca_result A list created by [pca_bigmatrix()] or
#'   [pca_stream_bigmatrix()] containing standard deviation and explained variance
#'   elements.
#' @param max_components Maximum number of components to display. Defaults to
#'   25 or the available number of components, whichever is smaller.
#' @param cumulative Logical flag indicating whether to overlay the cumulative
#'   explained variance line.
#' @param draw Logical; set to `FALSE` to return the prepared data without
#'   drawing a plot (useful for testing).
#' @param ... Additional parameters passed to [plot()].
#'
#' @return A list with `component`, `explained`, and `cumulative` vectors is
#'   returned invisibly. When `draw = TRUE`, the function produces a scree plot
#'   using base graphics.
#' @export
pca_plot_scree <- function(pca_result, max_components = 25L, cumulative = TRUE, draw = TRUE, ...) {
    if (is.null(pca_result$sdev) || is.null(pca_result$explained_variance)) {
        stop("`pca_result` must contain `sdev` and `explained_variance` components")
    }
    explained <- as.numeric(pca_result$explained_variance)
    total_components <- length(explained)
    if (total_components == 0) {
        stop("`pca_result$explained_variance` cannot be empty")
    }
    max_components <- as.integer(max_components)
    if (is.na(max_components) || max_components <= 0L) {
        stop("`max_components` must be a positive integer")
    }
    keep <- seq_len(min(total_components, max_components))
    explained_keep <- explained[keep]
    cumulative_keep <- if (!is.null(pca_result$cumulative_variance)) {
        as.numeric(pca_result$cumulative_variance)[keep]
    } else {
        cumsum(explained)[keep]
    }
    if (draw) {
        old_par <- graphics::par(no.readonly = TRUE)
        on.exit(graphics::par(old_par))
        graphics::plot(keep, explained_keep, type = "b", xlab = "Component",
             ylab = "Explained variance", ylim = c(0, max(1, max(explained_keep, cumulative_keep))),
             ...)
        if (cumulative) {
            graphics::lines(keep, cumulative_keep, type = "o", pch = 1, lty = 2, col = "grey40")
            graphics::legend("topright", legend = c("Explained", "Cumulative"),
                   lty = c(1, 2), pch = c(19, 1), col = c("black", "grey40"), bty = "n")
        }
    }
    invisible(list(component = keep, explained = explained_keep, cumulative = cumulative_keep))
}

extract_rows <- function(source, indices) {
    if (inherits(source, "big.matrix")) {
        return(source[indices, , drop = FALSE])
    }
    if (is.matrix(source) || is.data.frame(source)) {
        return(as.matrix(source[indices, , drop = FALSE]))
    }
    stop("Unsupported source type for row extraction")
}

#' Plot sampled PCA scores
#'
#' @description
#' Streams a subset of observations through the PCA rotation and plots their
#' scores on the requested components. Sampling keeps the drawn subset small so
#' graphics remain interpretable even when the source big matrix contains
#' millions of rows.
#'
#' @param x Either a `bigmemory::big.matrix`, a standard matrix, or a data
#'   frame.
#' @param rotation A rotation matrix such as `pca_result$rotation`.
#' @param center Optional centering vector. Use `numeric()` when no centering
#'   was applied.
#' @param scale Optional scaling vector. Use `numeric()` when no scaling was
#'   applied.
#' @param components Length-two integer vector selecting the principal
#'   components to display.
#' @param max_points Maximum number of observations to sample for the plot.
#' @param sample Strategy for selecting rows. `"uniform"` draws a random sample
#'   without replacement, whereas `"head"` takes the first `max_points` rows.
#' @param seed Optional seed to make the sampling reproducible.
#' @param draw Logical; set to `FALSE` to skip plotting and only return the
#'   sampled scores.
#' @param ... Additional graphical parameters forwarded to [plot()].
#'
#' @return A list containing `indices` (the sampled row indices) and `scores`
#'   (the corresponding score matrix) is returned invisibly. When `draw = TRUE`
#'   a scatter plot is produced.
#' @export
pca_plot_scores <- function(x, rotation, center = numeric(), scale = numeric(),
                            components = c(1L, 2L), max_points = 5000L,
                            sample = c("uniform", "head"), seed = NULL,
                            draw = TRUE, ...) {
    if (length(components) != 2L) {
        stop("`components` must be a vector of length two")
    }
    sample <- match.arg(sample)
    max_points <- as.integer(max_points)
    if (is.na(max_points) || max_points <= 0L) {
        stop("`max_points` must be a positive integer")
    }
    if (inherits(x, "big.matrix")) {
        total_rows <- nrow(x)
    } else if (is.matrix(x) || is.data.frame(x)) {
        total_rows <- nrow(x)
        x <- as.matrix(x)
    } else {
        stop("`x` must be a big.matrix, matrix, or data.frame")
    }
    if (total_rows == 0L) {
        stop("`x` must contain at least one row")
    }
    keep_n <- min(total_rows, max_points)
    if (!is.null(seed)) {
        set.seed(seed)
    }
    indices <- if (sample == "uniform" && keep_n < total_rows) {
        sort(sample.int(total_rows, keep_n, replace = FALSE))
    } else {
        seq_len(keep_n)
    }
    block <- extract_rows(x, indices)
    if (length(center)) {
        block <- sweep(block, 2, center, "-")
    }
    if (length(scale)) {
        safe_scale <- ifelse(scale == 0, 1, scale)
        block <- sweep(block, 2, safe_scale, "/")
    }
    rot <- as.matrix(rotation)
    comps <- as.integer(components)
    if (any(comps < 1) || any(comps > ncol(rot))) {
        stop("Requested components exceed available rotation columns")
    }
    scores <- block %*% rot[, comps, drop = FALSE]
    if (draw) {
        graphics::plot(scores[, 1], scores[, 2], xlab = sprintf("PC%d", comps[1]),
             ylab = sprintf("PC%d", comps[2]), ...)
    }
    invisible(list(indices = indices, scores = scores))
}

#' PCA biplot helper
#'
#' @description
#' Combines principal component scores and variable loadings in a single
#' scatter plot. The helper accepts both standard matrices and
#' `bigmemory::big.matrix` inputs, extracting only the requested component
#' columns. When `draw = TRUE`, the function scales the loadings to match the
#' score ranges, draws optional axes, overlays loading arrows, and labels
#' observations when requested.
#'
#' @param scores Matrix or `bigmemory::big.matrix` containing principal
#'   component scores with observations in rows and components in columns.
#' @param loadings Matrix or `bigmemory::big.matrix` of variable loadings whose
#'   columns correspond to principal components.
#' @param components Integer vector of length two selecting the components to
#'   display.
#' @param draw Logical; set to `FALSE` to return the prepared data without
#'   plotting.
#' @param draw_axes Logical; when `TRUE`, horizontal and vertical axes are
#'   drawn through the origin.
#' @param draw_arrows Logical; when `TRUE`, loading arrows are rendered.
#' @param label_points Logical; when `TRUE`, point labels derived from row names
#'   are drawn next to the scores.
#' @param ... Additional graphical parameters passed to [graphics::plot()].
#'
#' @return A list containing the selected `components`, extracted `scores`,
#'   original `loadings`, scaled loadings (`loadings_scaled`), and the applied
#'   `scale_factor`. The list is returned invisibly. When `draw = TRUE`, a biplot
#'   is produced using base graphics.
#' @export
pca_plot_biplot <- function(scores, loadings, components = c(1L, 2L), draw = TRUE,
                            draw_axes = TRUE, draw_arrows = TRUE,
                            label_points = FALSE, ...) {
    if (length(components) != 2L) {
        stop("`components` must contain exactly two indices")
    }
    comps <- as.integer(components)
    if (any(is.na(comps))) {
        stop("`components` must be coercible to integers")
    }

    original_comps <- comps

    ncol_scores <- ncol(scores)
    score_indices <- comps
    if (ncol_scores < max(comps)) {
        if (ncol_scores == length(comps)) {
            score_indices <- seq_len(ncol_scores)
        } else {
            stop("Requested components exceed the available score columns")
        }
    }

    ncol_loadings <- ncol(loadings)
    loading_indices <- comps
    if (ncol_loadings < max(comps)) {
        if (ncol_loadings == length(comps)) {
            loading_indices <- seq_len(ncol_loadings)
        } else {
            stop("Requested components exceed the available loading columns")
        }
    }

    extract_columns <- function(source, columns) {
        if (inherits(source, "big.matrix")) {
            return(as.matrix(source[, columns, drop = FALSE]))
        }
        if (is.matrix(source) || is.data.frame(source)) {
            mat <- as.matrix(source)
            return(mat[, columns, drop = FALSE])
        }
        stop("Inputs must be matrices, data.frames, or big.matrix objects")
    }

    score_mat <- extract_columns(scores, score_indices)
    loading_mat <- extract_columns(loadings, loading_indices)

    if (!nrow(score_mat) || !nrow(loading_mat)) {
        stop("`scores` and `loadings` must contain at least one row")
    }

    score_range <- apply(abs(score_mat), 2L, max, na.rm = TRUE)
    loading_range <- apply(abs(loading_mat), 2L, max, na.rm = TRUE)
    valid <- loading_range > 0 & is.finite(loading_range)
    if (!any(valid)) {
        scale_factor <- 1
    } else {
        ratios <- score_range[valid] / loading_range[valid]
        ratios[!is.finite(ratios) | ratios <= 0] <- NA_real_
        ratios <- ratios[!is.na(ratios)]
        scale_factor <- if (length(ratios)) min(ratios) else 1
    }
    if (!is.finite(scale_factor) || scale_factor <= 0) {
        scale_factor <- 1
    }
    loading_scaled <- loading_mat * scale_factor

    axis_labels <- sprintf("PC%d", original_comps)
    rownames_scores <- rownames(score_mat)
    if (is.null(rownames_scores)) {
        rownames_scores <- as.character(seq_len(nrow(score_mat)))
    }
    rownames_loadings <- rownames(loading_mat)
    if (is.null(rownames_loadings)) {
        rownames_loadings <- as.character(seq_len(nrow(loading_mat)))
    }

    if (draw) {
        combined <- rbind(score_mat, loading_scaled)
        x_range <- range(c(0, combined[, 1L]), finite = TRUE)
        y_range <- range(c(0, combined[, 2L]), finite = TRUE)
        if (!length(x_range) || !is.finite(x_range[1L])) {
            x_range <- c(-1, 1)
        }
        if (!length(y_range) || !is.finite(y_range[1L])) {
            y_range <- c(-1, 1)
        }
        graphics::plot(NA_real_, NA_real_, type = "n", xlab = axis_labels[1L],
            ylab = axis_labels[2L], xlim = x_range, ylim = y_range, ...)
        if (draw_axes) {
            graphics::abline(h = 0, v = 0, col = "grey80", lty = "dashed")
        }
        graphics::points(score_mat[, 1L], score_mat[, 2L], pch = 19)
        if (label_points) {
            graphics::text(score_mat[, 1L], score_mat[, 2L],
                labels = rownames_scores, pos = 3, cex = 0.7)
        }
        if (draw_arrows) {
            graphics::arrows(0, 0, loading_scaled[, 1L], loading_scaled[, 2L],
                length = 0.08, col = "red")
            graphics::text(loading_scaled[, 1L], loading_scaled[, 2L],
                labels = rownames_loadings, col = "red", pos = 4, cex = 0.8)
        }
    }

    invisible(list(
        components = original_comps,
        scores = score_mat,
        loadings = loading_mat,
        loadings_scaled = loading_scaled,
        scale_factor = scale_factor
    ))
}

#' Plot variable contributions
#' 
#' @description
#' Highlights the variables that contribute most to a selected principal
#' component. The helper works with dense matrices returned by
#' [pca_variable_contributions()] as well as with `bigmemory::big.matrix`
#' objects via sampling.
#'
#' @param contributions Contribution matrix where rows correspond to variables
#'   and columns to components.
#' @param component Integer index of the component to visualise.
#' @param top_n Number of variables with the largest absolute contribution to
#'   include in the bar plot.
#' @param draw Logical; set to `FALSE` to skip plotting.
#' @param ... Additional arguments passed to [barplot()].
#'
#' @return A data frame with the variables and their contributions is returned
#'   invisibly. When `draw = TRUE`, a bar plot of the top variables is produced.
#' @export
pca_plot_contributions <- function(contributions, component = 1L, top_n = 20L,
                                   draw = TRUE, ...) {
    if (inherits(contributions, "big.matrix")) {
        mat <- contributions[, , drop = FALSE]
    } else {
        mat <- as.matrix(contributions)
    }
    if (!ncol(mat)) {
        stop("`contributions` must have at least one column")
    }
    component <- as.integer(component)
    if (is.na(component) || component < 1L || component > ncol(mat)) {
        stop("`component` is out of range")
    }
    values <- mat[, component]
    ord <- order(abs(values), decreasing = TRUE)
    top_n <- as.integer(min(length(values), max(1L, top_n)))
    keep <- ord[seq_len(top_n)]
    variable_names <- rownames(mat)
    labels <- if (!is.null(variable_names)) variable_names[keep] else as.character(keep)
    data <- data.frame(
        variable = factor(labels, levels = labels),
        contribution = values[keep],
        stringsAsFactors = FALSE
    )
    if (draw) {
        graphics::barplot(height = data$contribution, names.arg = as.character(data$variable),
                las = 2, ylab = sprintf("Contribution to PC%d", component), ...)
    }
    invisible(data)
}

#' Plot a PCA correlation circle
#'
#' @description
#' Visualises the correlation between each variable and a pair of principal
#' components. The variables are projected onto the unit circle, where points
#' near the perimeter indicate strong correlation with the selected components.
#'
#' @param correlations Matrix or `bigmemory::big.matrix` containing variable
#'   correlations, typically produced by [pca_variable_correlations()].
#' @param components Length-two integer vector specifying the principal
#'   components to display.
#' @param labels Optional character vector specifying the labels to display for
#'   each variable. When `NULL`, the row names of `correlations` are used when
#'   available.
#' @param draw Logical; set to `FALSE` to return the prepared coordinates
#'   without plotting.
#' @param ... Additional graphical parameters passed to [graphics::plot()].
#'
#' @return A data frame with `variable`, `PCx`, and `PCy` columns representing the
#'   projected correlations, where `PCx`/`PCy` correspond to the requested
#'   component indices. The data frame is returned invisibly.
#' @export
pca_plot_correlation_circle <- function(correlations, components = c(1L, 2L),
                                        labels = NULL, draw = TRUE, ...) {
    if (length(components) != 2L) {
        stop("`components` must contain exactly two component indices")
    }
    comps <- as.integer(components)
    cor_mat <- if (inherits(correlations, "big.matrix")) {
        correlations[, , drop = FALSE]
    } else {
        as.matrix(correlations)
    }
    if (!ncol(cor_mat)) {
        stop("`correlations` must contain at least one component")
    }
    if (any(comps < 1L) || any(comps > ncol(cor_mat))) {
        stop("Requested components exceed the available correlation columns")
    }
    coords <- cor_mat[, comps, drop = FALSE]
    vars <- if (!is.null(labels)) {
        if (length(labels) != nrow(coords)) {
            stop("`labels` must have the same length as the number of variables")
        }
        labels
    } else if (!is.null(rownames(coords))) {
        rownames(coords)
    } else {
        seq_len(nrow(coords))
    }
    component_names <- paste0("PC", comps)
    coord_df <- as.data.frame(coords[, seq_along(comps), drop = FALSE])
    names(coord_df) <- component_names
    result <- cbind(
        data.frame(variable = vars, stringsAsFactors = FALSE),
        coord_df
    )
    x_col <- component_names[1]
    y_col <- component_names[2]

    if (draw) {
        graphics::plot(0, 0,
            type = "n",
            xlim = c(-1, 1),
            ylim = c(-1, 1),
            xlab = sprintf("Correlation with PC%d", comps[1]),
            ylab = sprintf("Correlation with PC%d", comps[2]),
            asp = 1,
            ...
        )
        graphics::symbols(0, 0, circles = 1, inches = FALSE, add = TRUE, fg = "grey70")
        graphics::abline(h = 0, v = 0, lty = 3, col = "grey80")
        graphics::arrows(0, 0, result[[x_col]], result[[y_col]], length = 0.08, col = "steelblue")
        graphics::text(result[[x_col]], result[[y_col]], labels = result$variable, pos = 4, cex = 0.8)
    }
    invisible(result)
}

# #' PCA biplot of scores and loadings
# #'
# #' @description
# #' Draws a classical PCA biplot that overlays observation scores and variable
# #' loadings on a shared coordinate system. Loadings are rescaled to fit the
# #' spread of the scores, mimicking the behaviour of [stats::biplot()].
# #'
# #' @param scores Matrix or `bigmemory::big.matrix` containing observation scores
# #'   on the principal components.
# #' @param loadings Matrix or `bigmemory::big.matrix` containing variable loadings
# #'   (typically the output of [pca_variable_loadings()]).
# #' @param components Length-two integer vector selecting the principal
# #'   components to display.
# #' @param scale Numeric scaling factor applied to the loading vectors after the
# #'   automatic rescaling step. Values larger than one emphasise loadings relative
# #'   to scores.
# #' @param draw Logical; set to `FALSE` to return the prepared coordinates without
# #'   drawing the plot.
# #' @param ... Additional graphical parameters passed to [graphics::plot()].
# #'
# #' @return A list containing the selected `scores`, original `loadings`,
# #'   rescaled `loadings_scaled`, and the applied `scale_factor`. The list is
# #'   returned invisibly.
# #' @export
#' pca_plot_biplot <- function(scores, loadings, components = c(1L, 2L), scale = 1,
#'                             draw = TRUE, ...) {
#'     if (length(components) != 2L) {
#'         stop("`components` must contain exactly two component indices")
#'     }
#'     comps <- as.integer(components)
#'     score_mat <- if (inherits(scores, "big.matrix")) {
#'         scores[, , drop = FALSE]
#'     } else {
#'         as.matrix(scores)
#'     }
#'     loading_mat <- if (inherits(loadings, "big.matrix")) {
#'         loadings[, , drop = FALSE]
#'     } else {
#'         as.matrix(loadings)
#'     }
#'     if (any(comps < 1L) || any(comps > ncol(score_mat))) {
#'         stop("`scores` do not contain the requested components")
#'     }
#'     if (any(comps < 1L) || any(comps > ncol(loading_mat))) {
#'         stop("`loadings` do not contain the requested components")
#'     }
#'     score_coords <- score_mat[, comps, drop = FALSE]
#'     loading_coords <- loading_mat[, comps, drop = FALSE]
#' 
#'     score_range <- apply(abs(score_coords), 2, function(x) if (length(x)) max(x, na.rm = TRUE) else 0)
#'     loading_range <- apply(abs(loading_coords), 2, function(x) if (length(x)) max(x, na.rm = TRUE) else 0)
#'     score_range[score_range == 0] <- 1
#'     loading_range[loading_range == 0] <- 1
#'     base_scale <- min(score_range / loading_range)
#'     if (!is.finite(base_scale) || base_scale <= 0) {
#'         base_scale <- 1
#'     }
#'     total_scale <- base_scale * scale
#'     loading_scaled <- loading_coords * total_scale
#' 
#'     result <- list(
#'         scores = score_coords,
#'         loadings = loading_coords,
#'         loadings_scaled = loading_scaled,
#'         scale_factor = total_scale
#'     )
#' 
#'     if (draw) {
#'         graphics::plot(score_coords[, 1], score_coords[, 2],
#'             xlab = sprintf("PC%d", comps[1]),
#'             ylab = sprintf("PC%d", comps[2]),
#'             ...
#'         )
#'         graphics::abline(h = 0, v = 0, lty = 3, col = "grey80")
#'         graphics::arrows(0, 0, loading_scaled[, 1], loading_scaled[, 2], length = 0.08, col = "firebrick")
#'         loading_labels <- rownames(loading_coords)
#'         if (is.null(loading_labels)) {
#'             loading_labels <- seq_len(nrow(loading_coords))
#'         }
#'         graphics::text(loading_scaled[, 1], loading_scaled[, 2], labels = loading_labels, col = "firebrick", pos = 4, cex = 0.8)
#'     }
#' 
#'     invisible(result)
#' }
