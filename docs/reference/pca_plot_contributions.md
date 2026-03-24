# Plot variable contributions

Highlights the variables that contribute most to a selected principal
component. The helper works with dense matrices returned by
[`pca_variable_contributions()`](https://fbertran.github.io/bigPCAcpp/reference/pca_bigmatrix.md)
as well as with
[`bigmemory::big.matrix`](https://rdrr.io/pkg/bigmemory/man/big.matrix.html)
objects via sampling.

## Usage

``` r
pca_plot_contributions(
  contributions,
  component = 1L,
  top_n = 20L,
  draw = TRUE,
  ...
)
```

## Arguments

- contributions:

  Contribution matrix where rows correspond to variables and columns to
  components.

- component:

  Integer index of the component to visualise.

- top_n:

  Number of variables with the largest absolute contribution to include
  in the bar plot.

- draw:

  Logical; set to `FALSE` to skip plotting.

- ...:

  Additional arguments passed to
  [`barplot()`](https://rdrr.io/r/graphics/barplot.html).

## Value

A data frame with the variables and their contributions is returned
invisibly. When `draw = TRUE`, a bar plot of the top variables is
produced.
