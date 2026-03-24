# Internal constructors and S3 methods for bigPCAcpp results

These helpers provide a light-weight S3 layer around the PCA outputs so
users can interact with them through familiar generics such as
[`summary()`](https://rdrr.io/r/base/summary.html) and
[`plot()`](https://rdrr.io/r/graphics/plot.default.html).

## Usage

``` r
new_bigpca_result(result, backend)
```
