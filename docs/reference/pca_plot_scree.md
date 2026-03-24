# Scree plot for principal component importance

Displays the proportion of variance explained by the leading principal
components. The function caps the number of displayed components to keep
the visualization legible on very high-dimensional problems.

## Usage

``` r
pca_plot_scree(
  pca_result,
  max_components = 25L,
  cumulative = TRUE,
  draw = TRUE,
  ...
)
```

## Arguments

- pca_result:

  A list created by
  [`pca_bigmatrix()`](https://fbertran.github.io/bigPCAcpp/reference/pca_bigmatrix.md)
  or
  [`pca_stream_bigmatrix()`](https://fbertran.github.io/bigPCAcpp/reference/pca_stream_bigmatrix.md)
  containing standard deviation and explained variance elements.

- max_components:

  Maximum number of components to display. Defaults to 25 or the
  available number of components, whichever is smaller.

- cumulative:

  Logical flag indicating whether to overlay the cumulative explained
  variance line.

- draw:

  Logical; set to `FALSE` to return the prepared data without drawing a
  plot (useful for testing).

- ...:

  Additional parameters passed to
  [`plot()`](https://rdrr.io/r/graphics/plot.default.html).

## Value

A list with `component`, `explained`, and `cumulative` vectors is
returned invisibly. When `draw = TRUE`, the function produces a scree
plot using base graphics.
