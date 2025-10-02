skip_if_not_installed("bigmemory")

create_big_matrix <- function(mat) {
  bm <- bigmemory::big.matrix(
    nrow = nrow(mat),
    ncol = ncol(mat),
    type = "double"
  )
  bm[,] <- mat
  bm
}

maybe_vec <- function(x) {
  if (is.null(x)) numeric() else x
}

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

test_that("svd_bigmatrix reproduces base svd results", {
  set.seed(2025)
  dense <- matrix(rnorm(24), nrow = 6, ncol = 4)
  bm <- create_big_matrix(dense)

  base <- svd(dense)
  res <- svd_bigmatrix(bm)

  expect_equal(res$d, base$d, tolerance = 1e-6)
  u_aligned <- align_columns(base$u, res$u)
  v_aligned <- align_columns(base$v, res$v)
  expect_equal(u_aligned, base$u, tolerance = 1e-6)
  expect_equal(v_aligned, base$v, tolerance = 1e-6)

  old_opts <- options(bigmemory.typecast.warning = FALSE)
  on.exit(options(old_opts), add = TRUE)

  float_bm <- bigmemory::big.matrix(nrow = nrow(dense), ncol = ncol(dense), type = "float")
  float_bm[,] <- dense
  float_res <- svd_bigmatrix(float_bm)
  u_float <- align_columns(base$u, float_res$u)
  v_float <- align_columns(base$v, float_res$v)
  expect_equal(float_res$d, base$d, tolerance = 1e-6)
  expect_equal(u_float, base$u, tolerance = 1e-6)
  expect_equal(v_float, base$v, tolerance = 1e-6)
})

test_that("svd_bigmatrix handles truncated outputs", {
  set.seed(404)
  dense <- matrix(rnorm(200), nrow = 20, ncol = 10)
  bm <- create_big_matrix(dense)

  res <- svd_bigmatrix(bm, nu = 0, nv = 0)
  base <- svd(dense, nu = 0, nv = 0)

  expect_equal(dim(res$u), c(nrow(dense), 0L))
  expect_equal(dim(res$v), c(ncol(dense), 0L))
  expect_equal(res$d, base$d, tolerance = 1e-6)

  pca_res <- pca_bigmatrix(bm, center = TRUE, scale = FALSE, ncomp = 5)
  expect_silent(pca_plot_scores(
    bm,
    pca_res$rotation,
    maybe_vec(pca_res$center),
    maybe_vec(pca_res$scale),
    max_points = 15,
    draw = FALSE
  ))
})

test_that("svd_bigmatrix rejects singular vector requests above the rank", {
  dense <- matrix(rnorm(60), nrow = 12, ncol = 5)
  bm <- create_big_matrix(dense)

  expect_error(
    svd_bigmatrix(bm, nu = min(nrow(dense), ncol(dense)) + 1, nv = 0),
    "`nu` cannot exceed min(nrow, ncol)",
    fixed = TRUE
  )
  expect_error(
    svd_bigmatrix(bm, nu = 0, nv = min(nrow(dense), ncol(dense)) + 1),
    "`nv` cannot exceed min(nrow, ncol)",
    fixed = TRUE
  )
})

test_that("pca_bigmatrix reproduces prcomp results", {
  set.seed(123)
  dense <- matrix(rnorm(30), nrow = 10, ncol = 3)
  bm <- create_big_matrix(dense)

  res <- pca_bigmatrix(bm, center = TRUE, scale = FALSE)
  expect_s3_class(res, c("bigpca_bigmemory", "bigpca"))
  pr <- prcomp(dense, center = TRUE, scale. = FALSE)

  expect_equal(res$nobs, nrow(dense))
  expect_null(res$scale)
  expect_equal(as.numeric(res$center), pr$center, tolerance = 1e-6)
  expect_equal(res$column_sd, apply(dense, 2, sd), tolerance = 1e-6)
  expect_equal(res$sdev, pr$sdev, tolerance = 1e-6)
  rotation_mat <- align_columns(unclass(pr$rotation), res$rotation[])
  expect_equal(rotation_mat, unclass(pr$rotation), tolerance = 1e-6)
  expect_equal(res$covariance, cov(dense), tolerance = 1e-6)

  explained <- (pr$sdev ^ 2) / sum(pr$sdev ^ 2)
  expect_equal(res$explained_variance, explained, tolerance = 1e-6)
  expect_equal(res$cumulative_variance, cumsum(explained), tolerance = 1e-6)

  summ <- summary(res)
  expect_s3_class(summ, "summary.bigpca")
  expect_equal(unname(summ$importance["Proportion of Variance", ]),
               as.numeric(res$explained_variance), tolerance = 1e-6)
  expect_silent(plot(res, type = "scree", max_components = 3, draw = FALSE))
  expect_silent(plot(res, type = "contributions", component = 1, top_n = 2, draw = FALSE))
  expect_silent(plot(res, type = "correlation_circle", draw = FALSE))

  scores <- pca_scores_bigmatrix(
    bm,
    res$rotation,
    maybe_vec(res$center),
    maybe_vec(res$scale)
  )
  aligned_scores <- align_columns(unname(pr$x), scores[])
  expect_equal(aligned_scores, unname(pr$x), tolerance = 1e-6)

  loadings <- pca_variable_loadings(res$rotation, res$sdev)
  expect_equal(loadings, sweep(res$rotation, 2, res$sdev, `*`), tolerance = 1e-6)

  correlations <- pca_variable_correlations(res$rotation, res$sdev, res$column_sd, res$scale)
  expected_cor <- sweep(loadings, 1, res$column_sd, `/`)
  expect_equal(correlations, expected_cor, tolerance = 1e-6)

  circle_data <- pca_plot_correlation_circle(correlations, draw = FALSE)
  expect_equal(nrow(circle_data), nrow(correlations))
  expect_equal(colnames(circle_data), c("variable", "PC1", "PC2"))
  expect_true(all(circle_data$PC1^2 + circle_data$PC2^2 <= 1 + 1e-8))

  correlations_big <- create_big_matrix(as.matrix(correlations))
  circle_big <- pca_plot_correlation_circle(correlations_big, draw = FALSE)
  expect_equal(circle_big$PC1, circle_data$PC1)
  expect_equal(circle_big$PC2, circle_data$PC2)

  contributions <- pca_variable_contributions(loadings)
  expect_equal(colSums(contributions), rep(1, ncol(contributions)), tolerance = 1e-6)
})

test_that("pca_variable_correlations stay bounded for scaled inputs", {
  set.seed(902)
  dense <- matrix(rnorm(200), nrow = 40, ncol = 5)
  bm <- create_big_matrix(dense)

  res <- pca_bigmatrix(bm, center = TRUE, scale = TRUE)
  loadings <- pca_variable_loadings(res$rotation, res$sdev)
  correlations <- pca_variable_correlations(res$rotation, res$sdev, res$column_sd, res$scale)

  corr_range <- range(correlations)
  expect_gte(corr_range[1], -1 - 1e-8)
  expect_lte(corr_range[2], 1 + 1e-8)
  expect_equal(correlations, loadings, tolerance = 1e-6)
})

test_that("pca_plot_biplot prepares scaled overlays", {
  set.seed(99)
  dense <- matrix(rnorm(80), nrow = 16, ncol = 5)
  bm <- create_big_matrix(dense)

  res <- pca_bigmatrix(bm, center = TRUE, scale = TRUE)
  scores_dense <- pca_scores_bigmatrix(
    bm,
    res$rotation,
    maybe_vec(res$center),
    maybe_vec(res$scale)
  )
  loadings <- pca_variable_loadings(res$rotation, res$sdev)

  biplot_data <- pca_plot_biplot(scores_dense, loadings, components = c(2, 1), draw = FALSE)
  expect_equal(biplot_data$components, c(2L, 1L))
  expect_equal(dim(biplot_data$scores), c(nrow(scores_dense), 2L))
  expect_equal(dim(biplot_data$loadings), c(nrow(loadings), 2L))
  expect_equal(dim(biplot_data$loadings_scaled), c(nrow(loadings), 2L))
  expect_true(biplot_data$scale_factor > 0)

  scores_big <- create_big_matrix(scores_dense)
  loadings_big <- create_big_matrix(loadings)
  big_data <- pca_plot_biplot(scores_big, loadings_big, draw = FALSE)
  expect_equal(big_data$components, c(1L, 2L))
  expect_equal(big_data$scores, biplot_data$scores[, c(2, 1)])
  expect_equal(dim(big_data$loadings_scaled), dim(biplot_data$loadings_scaled))
})

test_that("plot.bigpca draws biplots for arbitrary component pairs", {
  set.seed(642)
  mat <- matrix(rnorm(360), nrow = 40, ncol = 9)
  bm <- create_big_matrix(mat)

  res <- pca_bigmatrix(bm, center = TRUE, scale = TRUE, ncomp = 6)

  biplot_res <- plot(res, type = "biplot", data = mat, components = c(3, 4), draw = FALSE)

  expect_equal(biplot_res$components, c(3L, 4L))
  expect_equal(dim(biplot_res$scores), c(nrow(mat), 2L))
})

test_that("individual diagnostics compute contributions and cos2", {
  set.seed(111)
  dense <- matrix(rnorm(90), nrow = 18, ncol = 5)
  bm <- create_big_matrix(dense)

  res <- pca_bigmatrix(bm, center = TRUE, scale = TRUE)

  scores <- pca_scores_bigmatrix(
    bm,
    res$rotation,
    maybe_vec(res$center),
    maybe_vec(res$scale)
  )
  score_mat <- scores[]

  ind_contrib <- pca_individual_contributions(score_mat, res$sdev)
  ind_cos2 <- pca_individual_cos2(score_mat)

  expect_equal(colSums(ind_contrib), rep(1, ncol(ind_contrib)), tolerance = 1e-6)
  expect_equal(rowSums(ind_cos2), rep(1, nrow(ind_cos2)), tolerance = 1e-6)

  correlations <- pca_variable_correlations(res$rotation, res$sdev, res$column_sd, res$scale)
  var_cos2 <- pca_variable_cos2(correlations)
  expect_equal(var_cos2, correlations^2, tolerance = 1e-12)
})

test_that("supplementary diagnostics return expected shapes and values", {
  set.seed(222)
  dense <- matrix(rnorm(16 * 4), nrow = 16, ncol = 4)
  bm <- create_big_matrix(dense)

  res <- pca_bigmatrix(bm, center = TRUE, scale = TRUE)
  scores <- pca_scores_bigmatrix(
    bm,
    res$rotation,
    maybe_vec(res$center),
    maybe_vec(res$scale)
  )
  score_mat <- scores[]

  sup_ind <- matrix(rnorm(20 * ncol(dense)), nrow = 20, ncol = ncol(dense))
  colnames(sup_ind) <- colnames(dense)
  rownames(sup_ind) <- sprintf("sup_%02d", seq_len(nrow(sup_ind)))

  sup_ind_metrics <- pca_supplementary_individuals(
    sup_ind,
    res$rotation,
    res$sdev,
    center = maybe_vec(res$center),
    scale = maybe_vec(res$scale)
  )
  sup_ind_weighted <- pca_supplementary_individuals(
    sup_ind,
    res$rotation,
    res$sdev,
    center = maybe_vec(res$center),
    scale = maybe_vec(res$scale),
    total_weight = 5
  )

  expect_equal(dim(sup_ind_metrics$scores), c(nrow(sup_ind), ncol(res$rotation)))
  expect_equal(dim(sup_ind_metrics$contributions), dim(sup_ind_metrics$scores))
  expect_equal(dim(sup_ind_metrics$cos2), dim(sup_ind_metrics$scores))
  expect_equal(colSums(sup_ind_metrics$contributions), rep(1, ncol(res$rotation)), tolerance = 1e-6)
  expect_equal(colSums(sup_ind_weighted$contributions), rep(5, ncol(res$rotation)), tolerance = 1e-6)
  expect_equal(unname(rowSums(sup_ind_metrics$cos2)), rep(1, nrow(sup_ind)), tolerance = 1e-6)

  sup_vars <- matrix(rnorm(nrow(dense) * 3), nrow = nrow(dense), ncol = 3)
  colnames(sup_vars) <- sprintf("sup_var_%d", seq_len(ncol(sup_vars)))

  sup_var_metrics <- pca_supplementary_variables(sup_vars, score_mat, res$sdev)

  expect_equal(dim(sup_var_metrics$loadings), c(ncol(sup_vars), ncol(res$rotation)))
  expect_equal(dim(sup_var_metrics$correlations), dim(sup_var_metrics$loadings))
  expect_equal(dim(sup_var_metrics$contributions), dim(sup_var_metrics$loadings))
  expect_equal(dim(sup_var_metrics$cos2), dim(sup_var_metrics$loadings))
  expect_equal(
    sup_var_metrics$correlations,
    stats::cor(sup_vars, score_mat),
    tolerance = 1e-6
  )
  expect_equal(sup_var_metrics$cos2, sup_var_metrics$correlations^2, tolerance = 1e-12)
})

test_that("plot.bigpca dispatches to the biplot helper", {
  set.seed(1234)
  dense <- matrix(rnorm(60), nrow = 12, ncol = 5)
  bm <- create_big_matrix(dense)

  res <- pca_bigmatrix(bm, center = TRUE, scale = TRUE)

  direct_scores <- pca_plot_scores(
    bm,
    res$rotation,
    maybe_vec(res$center),
    maybe_vec(res$scale),
    draw = FALSE
  )
  loadings <- pca_variable_loadings(res$rotation, res$sdev)
  direct <- pca_plot_biplot(direct_scores$scores, loadings, draw = FALSE)

  expect_error(
    method_out <- plot(res, type = "biplot", data = bm, draw = FALSE),
    NA
  )
  expect_equal(method_out$components, direct$components)
  expect_equal(method_out$scores, direct$scores)
  expect_equal(method_out$loadings, direct$loadings)
})

test_that("pca_scores helpers accept NULL center and scale", {
  set.seed(456)
  dense <- matrix(rnorm(24), nrow = 8, ncol = 3)
  bm <- create_big_matrix(dense)

  res <- pca_bigmatrix(bm, center = FALSE, scale = FALSE)

  expect_null(res$center)
  expect_null(res$scale)

  expect_error(
    scores <- pca_scores_bigmatrix(bm, res$rotation, res$center, res$scale),
    NA
  )
  expect_type(scores, "double")

  scores_dest <- bigmemory::big.matrix(nrow = nrow(dense), ncol = ncol(res$rotation), type = "double")
  expect_error(
    pca_scores_stream_bigmatrix(
      bm,
      scores_dest,
      res$rotation,
      res$center,
      res$scale
    ),
    NA
  )
  expect_equal(scores_dest[, ], as.matrix(scores), tolerance = 1e-6)
})

test_that("external pointer inputs remain supported", {
  set.seed(321)
  dense <- matrix(rnorm(20), nrow = 5)
  bm <- create_big_matrix(dense)

  expect_error(pca_bigmatrix(bm@address, center = TRUE, scale = TRUE), NA)
  expect_error(pca_scores_bigmatrix(
    bm@address,
    matrix(1, ncol = 1, nrow = ncol(dense)),
    numeric(),
    numeric(),
    ncomp = 1
  ), NA)
})


skip_if_not_installed("bigmemory")

test_that("streaming big.matrix writers stream PCA results", {
  set.seed(321)
  dense <- matrix(rnorm(24), nrow = 8, ncol = 3)
  bm <- create_big_matrix(dense)

  rotation_dest <- bigmemory::big.matrix(nrow = ncol(dense), ncol = ncol(dense), type = "double")
  res <- pca_stream_bigmatrix(bm, xpRotation = rotation_dest, center = TRUE, scale = TRUE)
  expect_s3_class(res, c("bigpca_stream_bigmatrix", "bigpca"))
  expect_true(inherits(res$rotation_stream_bigmatrix, "externalptr"))
  expect_identical(res$rotation_stream_bigmatrix, rotation_dest@address)

  expect_equal(rotation_dest[, ], as.matrix(res$rotation), tolerance = 1e-6)

  scores_ref <- pca_scores_bigmatrix(
    bm,
    res$rotation,
    maybe_vec(res$center),
    maybe_vec(res$scale)
  )
  scores_dest <- bigmemory::big.matrix(nrow = nrow(dense), ncol = ncol(res$rotation), type = "double")
  pca_scores_stream_bigmatrix(
    bm,
    scores_dest,
    res$rotation,
    maybe_vec(res$center),
    maybe_vec(res$scale)
  )
  expect_equal(scores_dest[, ], as.matrix(scores_ref), tolerance = 1e-6)

  loadings_ref <- pca_variable_loadings(res$rotation, res$sdev)
  loadings_dest <- bigmemory::big.matrix(nrow = ncol(dense), ncol = ncol(res$rotation), type = "double")
  pca_variable_loadings_stream_bigmatrix(rotation_dest, res$sdev, loadings_dest)
  expect_equal(loadings_dest[, ], loadings_ref, tolerance = 1e-6)

  correlations_ref <- pca_variable_correlations(res$rotation, res$sdev, res$column_sd, res$scale)
  correlations_dest <- bigmemory::big.matrix(nrow = ncol(dense), ncol = ncol(res$rotation), type = "double")
  pca_variable_correlations_stream_bigmatrix(
    rotation_dest,
    res$sdev,
    res$column_sd,
    res$scale,
    correlations_dest
  )
  expect_equal(correlations_dest[, ], correlations_ref, tolerance = 1e-6)

  contributions_ref <- pca_variable_contributions(loadings_ref)
  contributions_dest <- bigmemory::big.matrix(nrow = ncol(dense), ncol = ncol(res$rotation), type = "double")
  pca_variable_contributions_stream_bigmatrix(loadings_dest, contributions_dest)
  expect_equal(contributions_dest[, ], contributions_ref, tolerance = 1e-6)
  expect_silent(summary(res))
  expect_silent(plot(res, type = "scree", draw = FALSE))
})

test_that("plotting helpers sample big data results", {
  skip_if_not_installed("bigmemory")

  set.seed(99)
  dense <- matrix(rnorm(120), nrow = 40, ncol = 3)
  bm <- create_big_matrix(dense)

  res <- pca_bigmatrix(bm, center = TRUE, scale = TRUE)

  scree <- pca_plot_scree(res, max_components = 2, draw = FALSE)
  expect_length(scree$component, 2)

  scores_info <- pca_plot_scores(bm, res$rotation, maybe_vec(res$center), maybe_vec(res$scale),
                                 max_points = 5, seed = 123, draw = FALSE)
  expect_equal(nrow(scores_info$scores), 5)
  expect_true(all(scores_info$indices %in% seq_len(nrow(bm))))

  loadings <- pca_variable_loadings(res$rotation, res$sdev)
  contributions <- pca_variable_contributions(loadings)
  contrib_data <- pca_plot_contributions(contributions, component = 1, top_n = 2, draw = FALSE)
  expect_equal(nrow(contrib_data), 2)

  correlations <- pca_variable_correlations(res$rotation, res$sdev, res$column_sd, res$scale)
  circle_data <- pca_plot_correlation_circle(correlations, components = c(1L, 2L), draw = FALSE)
  expect_equal(nrow(circle_data), ncol(dense))
  expect_named(circle_data, c("variable", "PC1", "PC2"))

  scores_full <- pca_scores_bigmatrix(bm, res$rotation, maybe_vec(res$center), maybe_vec(res$scale))
  biplot_info <- pca_plot_biplot(scores_full, loadings, components = c(1L, 2L), draw = FALSE)
  expect_true(all(c("scores", "loadings", "loadings_scaled", "scale_factor") %in% names(biplot_info)))
  expect_equal(ncol(biplot_info$scores), 2)

  expect_silent(plot(res, type = "correlation_circle", components = c(1, 2), draw = FALSE))
  expect_silent(plot(res, type = "biplot", components = c(1, 2), data = bm, draw = FALSE))
})

test_that("pca_robust uses median and MAD for centering and scaling", {
  set.seed(2024)
  dense <- matrix(rnorm(60), nrow = 12, ncol = 5)
  dense[1, 1] <- 100

  res <- pca_robust(dense, center = TRUE, scale = FALSE, ncomp = 3)
  expect_s3_class(res, c("bigpca_robust", "bigpca"))
  expect_equal(attr(res, "backend"), "robust")
  expect_equal(res$center, apply(dense, 2, stats::median))
  expect_null(res$scale)
  expect_equal(dim(res$rotation), c(ncol(dense), 3))
  expect_equal(dim(res$scores), c(nrow(dense), 3))

  scaled_res <- pca_robust(dense, center = TRUE, scale = TRUE, ncomp = 2)
  centered <- sweep(dense, 2, apply(dense, 2, stats::median), "-")
  expect_equal(scaled_res$scale, apply(centered, 2, stats::mad))
  expect_equal(dim(scaled_res$rotation), c(ncol(dense), 2))

  expect_error(pca_robust(dense, center = FALSE, scale = TRUE),
               "requires `center = TRUE`")
})

test_that("robust SVD down-weights high leverage observations", {
  set.seed(123)
  base_data <- matrix(rnorm(200), nrow = 40, ncol = 5)
  clean_res <- pca_robust(base_data, center = TRUE, scale = TRUE, ncomp = 2)

  contaminated <- base_data
  contaminated[1, ] <- contaminated[1, ] + 75

  robust_res <- pca_robust(contaminated, center = TRUE, scale = TRUE, ncomp = 2)

  expect_length(robust_res$robust_weights, nrow(contaminated))
  expect_true(max(robust_res$robust_weights) <= 1)
  expect_lt(robust_res$robust_weights[1], 1)

  classical_clean <- stats::prcomp(base_data, center = TRUE, scale. = TRUE)
  classical_cont <- stats::prcomp(contaminated, center = TRUE, scale. = TRUE)
  classical_explained <- (classical_cont$sdev^2) / sum(classical_cont$sdev^2)
  clean_classical_explained <- (classical_clean$sdev^2) / sum(classical_clean$sdev^2)
  diff_classical <- sum(abs(classical_explained - clean_classical_explained))
  diff_robust <- sum(abs(robust_res$explained_variance - clean_res$explained_variance))
  expect_lt(diff_robust, diff_classical)
})

test_that("robust PCA variance uses raw weights", {
  set.seed(456)
  base_data <- matrix(rnorm(150), nrow = 30, ncol = 5)
  contaminated <- base_data
  contaminated[1, ] <- contaminated[1, ] + 100

  res <- pca_robust(contaminated, center = TRUE, scale = TRUE, ncomp = 3)
  expect_true(any(res$robust_weights < 1))

  weights <- res$robust_weights
  scores <- res$scores
  weighted_norms <- colSums(sweep(scores^2, 1, weights, `*`))
  total_weight <- sum(weights)
  expected_sdev <- if (total_weight > 1) {
    sqrt(weighted_norms / (total_weight - 1))
  } else {
    sqrt(weighted_norms)
  }

  expect_equal(res$sdev, expected_sdev)

  eigenvalues <- expected_sdev^2
  total_var <- sum(eigenvalues)
  expected_explained <- if (total_var > 0) eigenvalues / total_var else rep(0, length(eigenvalues))

  expect_equal(res$explained_variance, expected_explained)
  expect_equal(res$cumulative_variance, cumsum(expected_explained))
})

test_that("C++ robust SVD matches the R reference implementation", {
  set.seed(321)
  x <- matrix(rnorm(160), nrow = 40, ncol = 4)
  x[1, ] <- x[1, ] + 50

  r_res <- bigPCAcpp:::svd_robust_R(x, ncomp = 3, max_iter = 10)
  cpp_res <- bigPCAcpp:::svd_robust(x, ncomp = 3, max_iter = 10)

  expect_equal(cpp_res$d, r_res$d, tolerance = 1e-6)
  expect_equal(cpp_res$weights, r_res$weights, tolerance = 1e-6)
  expect_equal(cpp_res$iterations, r_res$iterations)

  u_aligned <- align_columns(r_res$u, cpp_res$u)
  v_aligned <- align_columns(r_res$v, cpp_res$v)
  expect_equal(u_aligned, r_res$u, tolerance = 1e-5)
  expect_equal(v_aligned, r_res$v, tolerance = 1e-5)
})
