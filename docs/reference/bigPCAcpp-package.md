# bigPCAcpp: Principal Component Analysis for bigmemory Matrices

The **bigPCAcpp** package provides high-performance principal component
analysis routines that work directly with
[`bigmemory::big.matrix`](https://rdrr.io/pkg/bigmemory/man/big.matrix.html)
objects. Data are streamed through BLAS and LAPACK kernels so large,
file-backed matrices can be analysed without materialising dense copies
in R. Companion helpers compute scores, loadings, correlations, and
contributions, including streaming variants that write results to
[`bigmemory::big.matrix`](https://rdrr.io/pkg/bigmemory/man/big.matrix.html)
destinations used by file-based pipelines.

## See also

[`pca_bigmatrix()`](https://fbertran.github.io/bigPCAcpp/reference/pca_bigmatrix.md),
[`pca_stream_bigmatrix()`](https://fbertran.github.io/bigPCAcpp/reference/pca_stream_bigmatrix.md)

## Author

**Maintainer**: Frederic Bertrand <frederic.bertrand@lecnam.net>

## Examples

``` r
# \donttest{
library(bigmemory)
mat <- as.big.matrix(matrix(rnorm(20), nrow = 5))
result <- pca_bigmatrix(mat)
result$sdev
#> [1] 2.0348820 1.0310207 0.5643723 0.1060603
# }
```
