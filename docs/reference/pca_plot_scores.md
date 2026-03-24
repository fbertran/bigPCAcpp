# Plot sampled PCA scores

Streams a subset of observations through the PCA rotation and plots
their scores on the requested components. Sampling keeps the drawn
subset small so graphics remain interpretable even when the source big
matrix contains millions of rows.

## Usage

``` r
pca_plot_scores(
  x,
  rotation,
  center = numeric(),
  scale = numeric(),
  components = c(1L, 2L),
  max_points = 5000L,
  sample = c("uniform", "head"),
  seed = NULL,
  draw = TRUE,
  ...
)
```

## Arguments

- x:

  Either a
  [`bigmemory::big.matrix`](https://rdrr.io/pkg/bigmemory/man/big.matrix.html),
  a standard matrix, or a data frame.

- rotation:

  A rotation matrix such as `pca_result$rotation`.

- center:

  Optional centering vector. Use
  [`numeric()`](https://rdrr.io/r/base/numeric.html) when no centering
  was applied.

- scale:

  Optional scaling vector. Use
  [`numeric()`](https://rdrr.io/r/base/numeric.html) when no scaling was
  applied.

- components:

  Length-two integer vector selecting the principal components to
  display.

- max_points:

  Maximum number of observations to sample for the plot.

- sample:

  Strategy for selecting rows. `"uniform"` draws a random sample without
  replacement, whereas `"head"` takes the first `max_points` rows.

- seed:

  Optional seed to make the sampling reproducible.

- draw:

  Logical; set to `FALSE` to skip plotting and only return the sampled
  scores.

- ...:

  Additional graphical parameters forwarded to
  [`plot()`](https://rdrr.io/r/graphics/plot.default.html).

## Value

A list containing `indices` (the sampled row indices) and `scores` (the
corresponding score matrix) is returned invisibly. When `draw = TRUE` a
scatter plot is produced.
