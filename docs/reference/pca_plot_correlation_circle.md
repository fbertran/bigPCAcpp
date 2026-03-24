# Plot a PCA correlation circle

Visualises the correlation between each variable and a pair of principal
components. The variables are projected onto the unit circle, where
points near the perimeter indicate strong correlation with the selected
components.

## Usage

``` r
pca_plot_correlation_circle(
  correlations,
  components = c(1L, 2L),
  labels = NULL,
  draw = TRUE,
  ...
)
```

## Arguments

- correlations:

  Matrix or
  [`bigmemory::big.matrix`](https://rdrr.io/pkg/bigmemory/man/big.matrix.html)
  containing variable correlations, typically produced by
  [`pca_variable_correlations()`](https://fbertran.github.io/bigPCAcpp/reference/pca_bigmatrix.md).

- components:

  Length-two integer vector specifying the principal components to
  display.

- labels:

  Optional character vector specifying the labels to display for each
  variable. When `NULL`, the row names of `correlations` are used when
  available.

- draw:

  Logical; set to `FALSE` to return the prepared coordinates without
  plotting.

- ...:

  Additional graphical parameters passed to
  [`graphics::plot()`](https://rdrr.io/r/graphics/plot.default.html).

## Value

A data frame with `variable`, `PCx`, and `PCy` columns representing the
projected correlations, where `PCx`/`PCy` correspond to the requested
component indices. The data frame is returned invisibly.
