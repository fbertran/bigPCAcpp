#' bigPCAcpp: Principal Component Analysis for bigmemory Matrices
#'
#' The **bigPCAcpp** package provides high-performance principal component analysis
#' routines that work directly with [`bigmemory::big.matrix`] objects. Data are
#' streamed through BLAS and LAPACK kernels so large, file-backed matrices can be
#' analysed without materialising dense copies in R. Companion helpers compute
#' scores, loadings, correlations, and contributions, including streaming
#' variants that write results to `bigmemory::big.matrix` destinations used by
#' file-based pipelines.
#'
#' @aliases bigPCAcpp-package bigPCAcpp NULL
#' 
#' @seealso [pca_bigmatrix()], [pca_stream_bigmatrix()]
#' 
#' @examples
#' \dontrun{
#' library(bigmemory)
#' mat <- as.big.matrix(matrix(rnorm(20), nrow = 5))
#' result <- pca_bigmatrix(mat)
#' result$sdev
#' }
"_PACKAGE"

#' @useDynLib bigPCAcpp, .registration = TRUE
#' @importFrom Rcpp evalCpp
#' @importFrom stats rnorm
NULL