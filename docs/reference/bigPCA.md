# BigPCA result objects

Results returned by
[`pca_bigmatrix()`](https://fbertran.github.io/bigPCAcpp/reference/pca_bigmatrix.md),
[`pca_stream_bigmatrix()`](https://fbertran.github.io/bigPCAcpp/reference/pca_stream_bigmatrix.md),
and
[`pca_robust()`](https://fbertran.github.io/bigPCAcpp/reference/pca_robust.md)
inherit from the `bigpca` class. The objects store the component
standard deviations, rotation/loadings, and optional scores while
recording which computational backend produced them. Standard S3
generics such as [`summary()`](https://rdrr.io/r/base/summary.html) and
[`plot()`](https://rdrr.io/r/graphics/plot.default.html) are implemented
for convenience.

`bigpca` objects are lists produced by
[`pca_bigmatrix()`](https://fbertran.github.io/bigPCAcpp/reference/pca_bigmatrix.md),
[`pca_stream_bigmatrix()`](https://fbertran.github.io/bigPCAcpp/reference/pca_stream_bigmatrix.md),
[`pca_robust()`](https://fbertran.github.io/bigPCAcpp/reference/pca_robust.md),
and related helpers. They mirror the structure of base R's
[`prcomp()`](https://rdrr.io/r/stats/prcomp.html) outputs while tracking
additional metadata for large-scale and streaming computations.

\#' @seealso
[`pca_bigmatrix()`](https://fbertran.github.io/bigPCAcpp/reference/pca_bigmatrix.md),
[`pca_stream_bigmatrix()`](https://fbertran.github.io/bigPCAcpp/reference/pca_stream_bigmatrix.md),
[`pca_robust()`](https://fbertran.github.io/bigPCAcpp/reference/pca_robust.md),
[`pca_plot_scree()`](https://fbertran.github.io/bigPCAcpp/reference/pca_plot_scree.md),
[`pca_plot_scores()`](https://fbertran.github.io/bigPCAcpp/reference/pca_plot_scores.md),
[`pca_plot_contributions()`](https://fbertran.github.io/bigPCAcpp/reference/pca_plot_contributions.md),
[`pca_plot_correlation_circle()`](https://fbertran.github.io/bigPCAcpp/reference/pca_plot_correlation_circle.md),
and
[`pca_plot_biplot()`](https://fbertran.github.io/bigPCAcpp/reference/pca_plot_biplot.md).

## Components

- `sdev`:

  Numeric vector of component standard deviations.

- `rotation`:

  Numeric matrix whose columns contain the variable loadings (principal
  axes).

- `center`, `scale`:

  Optional numeric vectors describing the centring and scaling applied
  to each variable when fitting the model.

- `scores`:

  Optional numeric matrix of principal component scores when computed
  alongside the decomposition.

- `column_sd`:

  Numeric vector of marginal standard deviations for each input
  variable.

- `eigenvalues`:

  Numeric vector of eigenvalues associated with the retained components.

- `explained_variance`, `cumulative_variance`:

  Numeric vectors summarising the fraction of variance explained by
  individual components and the corresponding cumulative totals.

- `covariance`:

  Sample covariance matrix used to derive the components.

- `nobs`:

  Number of observations used in the decomposition.

The class also records the computation backend via `attr(x, "backend")`,
enabling downstream methods to adjust their behaviour for streamed or
robust results.

## See also

[`pca_bigmatrix()`](https://fbertran.github.io/bigPCAcpp/reference/pca_bigmatrix.md),
[`pca_stream_bigmatrix()`](https://fbertran.github.io/bigPCAcpp/reference/pca_stream_bigmatrix.md),
[`summary.bigpca()`](https://fbertran.github.io/bigPCAcpp/reference/pca_bigmatrix.md),
[`print.summary.bigpca()`](https://fbertran.github.io/bigPCAcpp/reference/pca_bigmatrix.md),
[`plot.bigpca()`](https://fbertran.github.io/bigPCAcpp/reference/pca_bigmatrix.md)
