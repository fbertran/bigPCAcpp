suppressPackageStartupMessages({
  library(bigmemory)
  if (requireNamespace("bigPCAcpp", quietly = TRUE)) {
    library(bigPCAcpp)
  } else {
    if (!requireNamespace("pkgload", quietly = TRUE)) {
      stop("bigPCAcpp must be installed or pkgload must be available", call. = FALSE)
    }
    pkgload::load_all(".")
  }
})

if (!requireNamespace("bench", quietly = TRUE)) {
  stop("Package 'bench' is required for this benchmark", call. = FALSE)
}

if (!requireNamespace("irlba", quietly = TRUE)) {
  stop("Package 'irlba' is required for this benchmark", call. = FALSE)
}

set.seed(2025)
bm <- bigmemory::big.matrix(nrow = 2500L, ncol = 40L, type = "double")
m <- matrix(rnorm(2500L * 40L), nrow = 2500L)
bm[,] <- m

benchmark_irlba_results <- bench::mark(
  bigpca = pca_bigmatrix(bm, center = TRUE, scale = TRUE, ncomp = 4)$dev,
  irlba = irlba::prcomp_irlba(
    m,
    n = 4,
    center = TRUE,
    scale. = TRUE
  )$dev,
  min_iterations = 200L
)

print(benchmark_irlba_results)
