# Supplementary variable diagnostics

Compute loadings, correlations, contributions, and cos^2 values for
supplementary variables (columns) given component scores for the active
individuals.

## Usage

``` r
pca_supplementary_variables(data, scores, sdev, center = NULL)
```

## Arguments

- data:

  Matrix-like object whose columns correspond to supplementary variables
  measured on the active individuals.

- scores:

  Numeric matrix of component scores for the active individuals.

- sdev:

  Numeric vector of component standard deviations associated with
  `scores`.

- center:

  Optional numeric vector specifying the centring to apply to each
  supplementary variable. When `NULL`, column means of `data` are used.

## Value

A list with elements `loadings`, `correlations`, `contributions`, and
`cos2`.
