# Supplementary individual diagnostics

Compute principal component scores and quality metrics for supplementary
individuals (rows) projected into an existing PCA solution.

## Usage

``` r
pca_supplementary_individuals(
  data,
  rotation,
  sdev,
  center = NULL,
  scale = NULL,
  total_weight = NA_real_
)
```

## Arguments

- data:

  Matrix-like object whose rows correspond to supplementary individuals
  and columns to the original variables.

- rotation:

  Rotation matrix from the PCA model (e.g. the `rotation` element of a
  [`bigpca`](https://fbertran.github.io/bigPCAcpp/reference/bigpca.md)
  result).

- sdev:

  Numeric vector of component standard deviations associated with
  `rotation`.

- center:

  Optional numeric vector giving the centring applied to each variable
  when fitting the PCA. Defaults to zero centring.

- scale:

  Optional numeric vector describing the scaling applied to each
  variable when fitting the PCA. When `NULL`, no scaling is applied.

- total_weight:

  Optional positive scalar passed to
  [`pca_individual_contributions()`](https://fbertran.github.io/bigPCAcpp/reference/pca_bigmatrix.md)
  when computing contributions. When left as `NA` (the default), the
  resulting contributions for each component are normalised to sum to
  one across supplementary individuals. Supplying a value bypasses this
  normalisation and delegates the scaling to
  [`pca_individual_contributions()`](https://fbertran.github.io/bigPCAcpp/reference/pca_bigmatrix.md).

## Value

A list with elements `scores`, `contributions`, and `cos2`.
