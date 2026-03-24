# Plot PCA diagnostics for big data workflows

These helpers visualise the results returned by
[`pca_bigmatrix()`](https://fbertran.github.io/bigPCAcpp/reference/pca_bigmatrix.md)
and its companions without requiring users to materialise dense
intermediate structures. Each plotting function optionally samples the
inputs so the default output remains responsive even when the underlying
big matrix spans millions of observations.

## See also

[`pca_bigmatrix()`](https://fbertran.github.io/bigPCAcpp/reference/pca_bigmatrix.md),
[`pca_variable_loadings()`](https://fbertran.github.io/bigPCAcpp/reference/pca_bigmatrix.md),
[`pca_variable_contributions()`](https://fbertran.github.io/bigPCAcpp/reference/pca_bigmatrix.md)
