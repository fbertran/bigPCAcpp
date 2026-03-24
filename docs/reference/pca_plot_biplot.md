# PCA biplot helper

Combines principal component scores and variable loadings in a single
scatter plot. The helper accepts both standard matrices and
[`bigmemory::big.matrix`](https://rdrr.io/pkg/bigmemory/man/big.matrix.html)
inputs, extracting only the requested component columns. When
`draw = TRUE`, the function scales the loadings to match the score
ranges, draws optional axes, overlays loading arrows, and labels
observations when requested.

## Usage

``` r
pca_plot_biplot(
  scores,
  loadings,
  components = c(1L, 2L),
  draw = TRUE,
  draw_axes = TRUE,
  draw_arrows = TRUE,
  label_points = FALSE,
  ...
)
```

## Arguments

- scores:

  Matrix or
  [`bigmemory::big.matrix`](https://rdrr.io/pkg/bigmemory/man/big.matrix.html)
  containing principal component scores with observations in rows and
  components in columns.

- loadings:

  Matrix or
  [`bigmemory::big.matrix`](https://rdrr.io/pkg/bigmemory/man/big.matrix.html)
  of variable loadings whose columns correspond to principal components.

- components:

  Integer vector of length two selecting the components to display.

- draw:

  Logical; set to `FALSE` to return the prepared data without plotting.

- draw_axes:

  Logical; when `TRUE`, horizontal and vertical axes are drawn through
  the origin.

- draw_arrows:

  Logical; when `TRUE`, loading arrows are rendered.

- label_points:

  Logical; when `TRUE`, point labels derived from row names are drawn
  next to the scores.

- ...:

  Additional graphical parameters passed to
  [`graphics::plot()`](https://rdrr.io/r/graphics/plot.default.html).

## Value

A list containing the selected `components`, extracted `scores`,
original `loadings`, scaled loadings (`loadings_scaled`), and the
applied `scale_factor`. The list is returned invisibly. When
`draw = TRUE`, a biplot is produced using base graphics.
