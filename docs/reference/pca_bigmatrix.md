# Principal component analysis for `bigmemory::big.matrix` inputs

Perform principal component analysis (PCA) directly on a
[`bigmemory::big.matrix`](https://rdrr.io/pkg/bigmemory/man/big.matrix.html)
without copying the data into R memory. The exported helpers mirror the
structure of base R's [`prcomp()`](https://rdrr.io/r/stats/prcomp.html)
while avoiding the need to materialise large matrices.

## Usage

``` r
resolve_big_pointer(x, arg, allow_null = FALSE)

pca_scores_bigmatrix(
  xpMat,
  rotation,
  center,
  scale,
  ncomp = -1L,
  block_size = 1024L
)

pca_variable_loadings(rotation, sdev)

pca_variable_correlations(rotation, sdev, column_sd, scale = NULL)

pca_variable_contributions(loadings)

pca_individual_contributions(scores, sdev, total_weight = NA_real_)

pca_individual_cos2(scores)

pca_variable_cos2(correlations)

# S3 method for class 'bigpca'
summary(object, ...)

# S3 method for class 'summary.bigpca'
print(x, digits = max(3, getOption("digits") - 3), ...)

# S3 method for class 'bigpca'
plot(
  x,
  y,
  type = c("scree", "contributions", "correlation_circle", "biplot"),
  max_components = 25L,
  component = 1L,
  top_n = 20L,
  components = c(1L, 2L),
  data = NULL,
  draw = TRUE,
  ...
)
```

## Arguments

- x:

  A `summary.bigpca` object.

- arg:

  Character string naming the argument being validated. Used to
  construct informative error messages.

- allow_null:

  Logical flag indicating whether `NULL` is accepted for the argument.
  When `TRUE`, a `NULL` input is returned unchanged.

- xpMat:

  Either a
  [`bigmemory::big.matrix`](https://rdrr.io/pkg/bigmemory/man/big.matrix.html)
  or an external pointer such as `mat@address` that references the
  source `big.matrix`.

- rotation:

  A rotation matrix such as the `rotation` element returned by
  `pca_bigmatrix()`.

- center:

  For `pca_scores_bigmatrix()`, a numeric vector of column means
  (optional).

- scale:

  Optional numeric vector of scaling factors returned by
  `pca_bigmatrix()`. If supplied, it indicates the PCA was performed on
  standardised variables.

- ncomp:

  Number of components to retain. Use a non-positive value to keep all
  components returned by the decomposition.

- block_size:

  Number of rows to process per block when streaming data through BLAS
  kernels. Larger values improve throughput at the cost of additional
  memory.

- sdev:

  A numeric vector of component standard deviations, typically the
  `sdev` element from `pca_bigmatrix()`.

- column_sd:

  A numeric vector with the marginal standard deviation of each original
  variable. When `scale` is supplied, correlations are computed on the
  standardised scale without rescaling by `column_sd`.

- loadings:

  A numeric matrix such as the result of `pca_variable_loadings()`.

- scores:

  For `pca_individual_contributions()` and `pca_individual_cos2()`, a
  numeric matrix of component scores where rows correspond to
  observations and columns to components.

- total_weight:

  Optional positive scalar giving the effective number of observations
  to use when computing contributions. Defaults to the number of rows in
  `scores`.

- correlations:

  For `pca_variable_cos2()`, a numeric matrix of correlations between
  variables and components.

- object:

  A [`bigpca`](https://fbertran.github.io/bigPCAcpp/reference/bigpca.md)
  object created by `pca_bigmatrix()`,
  [`pca_stream_bigmatrix()`](https://fbertran.github.io/bigPCAcpp/reference/pca_stream_bigmatrix.md),
  or related helpers.

- ...:

  Additional arguments passed to plotting helpers.

- digits:

  Number of significant digits to display when printing importance
  metrics.

- y:

  Currently unused.

- type:

  The plot to draw. Options include "scree" (variance explained),
  "contributions" (top contributing variables), "correlation_circle"
  (variable correlations with selected components), and "biplot" (joint
  display of scores and loadings).

- max_components:

  Maximum number of components to display in scree plots.

- component:

  Component index to highlight when drawing contribution plots.

- top_n:

  Number of variables to display in contribution plots.

- components:

  Length-two integer vector selecting the components for correlation
  circle and biplot views.

- data:

  Optional data source (matrix, data frame,
  [`bigmemory::big.matrix`](https://rdrr.io/pkg/bigmemory/man/big.matrix.html),
  or external pointer) used to compute scores for biplots when
  `x$scores` is unavailable.

- draw:

  Logical; if `FALSE`, return the data prepared for the selected plot
  instead of drawing it.

## Value

For `pca_bigmatrix()`, a `bigpca` object mirroring a `prcomp` result
with elements `sdev`, `rotation`, optional `center` and `scale` vectors,
`column_sd`, `eigenvalues`, `explained_variance`, `cumulative_variance`,
and the sample covariance matrix. The object participates in S3 generics
such as [`summary()`](https://rdrr.io/r/base/summary.html) and
[`plot()`](https://rdrr.io/r/graphics/plot.default.html).

A numeric matrix of scores with rows corresponding to observations and
columns to retained components.

A numeric matrix containing variable loadings for each component.

A numeric matrix of correlations between variables and components.

A numeric matrix where each entry represents the contribution of a
variable to a component.

For `summary.bigpca()`, a `summary.bigpca` object containing component
importance measures.

## Functions

- `pca_scores_bigmatrix()`: Project observations into principal
  component space while streaming from a `big.matrix`.

- `pca_variable_loadings()`: Compute variable loadings (covariances
  between original variables and components).

- `pca_variable_correlations()`: Compute variable-component correlations
  given column standard deviations.

- `pca_variable_contributions()`: Derive the relative contribution of
  each variable to the retained components.

- `pca_individual_contributions()`: Compute the relative contribution of
  individual observations to each component.

- `pca_individual_cos2()`: Compute squared cosine values measuring the
  quality of representation for individual observations.

- `pca_variable_cos2()`: Compute squared cosine values measuring the
  quality of representation for variables.

- `summary(bigpca)`: Summarise the component importance metrics for a
  [`bigpca`](https://fbertran.github.io/bigPCAcpp/reference/bigpca.md)
  result.

- `print(summary.bigpca)`: Print the component importance summary
  produced by `summary.bigpca()`.

- `plot(bigpca)`: Visualise PCA diagnostics such as scree, correlation
  circle, contribution, and biplot displays.

## See also

[bigpca](https://fbertran.github.io/bigPCAcpp/reference/bigpca.md),
`pca_scores_bigmatrix()`, `pca_variable_loadings()`,
`pca_variable_correlations()`, `pca_variable_contributions()`, and the
streaming variants
[`pca_stream_bigmatrix()`](https://fbertran.github.io/bigPCAcpp/reference/pca_stream_bigmatrix.md)
and companions.

[bigpca](https://fbertran.github.io/bigPCAcpp/reference/bigpca.md)

[bigpca](https://fbertran.github.io/bigPCAcpp/reference/bigpca.md)

[bigpca](https://fbertran.github.io/bigPCAcpp/reference/bigpca.md)

## Examples

``` r
set.seed(123)
mat <- bigmemory::as.big.matrix(matrix(rnorm(40), nrow = 10))
pca <- pca_bigmatrix(mat, center = TRUE, scale = TRUE, ncomp = 3)
scores <- pca_scores_bigmatrix(mat, pca$rotation, pca$center, pca$scale, ncomp = 3)
loadings <- pca_variable_loadings(pca$rotation, pca$sdev)
correlations <- pca_variable_correlations(pca$rotation, pca$sdev, pca$column_sd, pca$scale)
contributions <- pca_variable_contributions(loadings)
list(scores = scores, loadings = loadings, correlations = correlations,
     contributions = contributions)
#> $scores
#>              [,1]        [,2]        [,3]
#>  [1,]  0.54325339 -0.91493729 -0.41315014
#>  [2,] -0.77428157 -0.84150780  0.40023111
#>  [3,]  1.77362948  0.79996096  0.05039649
#>  [4,]  0.61761772  0.59591822 -0.57031347
#>  [5,]  0.23218256  0.95503144 -0.72268699
#>  [6,]  2.66813521 -0.45268284  0.33775240
#>  [7,] -0.08448606  0.81203899  1.15087152
#>  [8,] -2.43891736  0.54722825 -0.96050148
#>  [9,] -0.40827305 -1.59689839 -0.36789931
#> [10,] -2.12886032  0.09584844  1.09529988
#> 
#> $loadings
#>            [,1]       [,2]       [,3]
#> [1,]  0.8567704  0.2711956  0.3439932
#> [2,]  0.7678169 -0.5050498  0.3324763
#> [3,] -0.7573969  0.3708264  0.5203357
#> [4,]  0.7581891  0.5754452 -0.2056261
#> 
#> $correlations
#>            [,1]       [,2]       [,3]
#> [1,]  0.8567704  0.2711956  0.3439932
#> [2,]  0.7678169 -0.5050498  0.3324763
#> [3,] -0.7573969  0.3708264  0.5203357
#> [4,]  0.7581891  0.5754452 -0.2056261
#> 
#> $contributions
#>           [,1]       [,2]       [,3]
#> [1,] 0.2969361 0.09224839 0.21836250
#> [2,] 0.2384786 0.31993524 0.20398568
#> [3,] 0.2320498 0.17247851 0.49962661
#> [4,] 0.2325354 0.41533786 0.07802521
#> 
```
