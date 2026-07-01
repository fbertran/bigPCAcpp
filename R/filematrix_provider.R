#' Filematrix Row-Block Provider
#'
#' `make_filematrix_row_provider()` creates a minimal block-access adapter for
#' `filematrix` objects. The provider is intended for block-wise sequential
#' access by benchmark and future streaming PCA backends. It avoids adding a
#' true filematrix PCA backend while exposing the same row-block contract used
#' by external benchmark layers.
#'
#' @param x,X Object to test or wrap.
#' @param chunk_size Suggested row chunk size for downstream consumers.
#'
#' @return `is_filematrix_object()` returns a logical scalar. The provider is a
#'   list with `nrow()`, `ncol()`, `get_rows()`, `get_cols()`, `get_block()`, and
#'   `storage_type = "filematrix"`.
#' @name filematrix_provider
#' @export
is_filematrix_object <- function(x) {
  inherits(x, "filematrix")
}

#' @rdname filematrix_provider
#' @export
make_filematrix_row_provider <- function(X, chunk_size = 1024L) {
  if (!requireNamespace("filematrix", quietly = TRUE)) {
    stop(
      "make_filematrix_row_provider() requires the optional filematrix package. ",
      "Install filematrix or use another storage backend.",
      call. = FALSE
    )
  }

  if (!is_filematrix_object(X)) {
    if (inherits(X, "sparseMatrix")) {
      stop("X must be a filematrix object; sparse matrix classes are not supported.", call. = FALSE)
    }
    stop("X must be a filematrix object.", call. = FALSE)
  }

  nr <- nrow(X)
  nc <- ncol(X)
  .validate_provider_dim(nr, "nrow(X)")
  .validate_provider_dim(nc, "ncol(X)")

  chunk_size <- .as_positive_scalar_integer(chunk_size, "chunk_size")
  .check_filematrix_numeric_storage(X)

  get_block <- function(i0, i1, j0, j1) {
    i0 <- .as_positive_scalar_integer(i0, "i0")
    i1 <- .as_positive_scalar_integer(i1, "i1")
    j0 <- .as_positive_scalar_integer(j0, "j0")
    j1 <- .as_positive_scalar_integer(j1, "j1")

    .validate_closed_interval(i0, i1, nr, "row")
    .validate_closed_interval(j0, j1, nc, "column")

    rows <- seq.int(i0, i1)
    cols <- seq.int(j0, j1)
    out <- X[rows, cols]
    .as_dense_double_matrix(out, nrow = length(rows), ncol = length(cols))
  }

  provider <- list(
    nrow = function() nr,
    ncol = function() nc,
    get_rows = function(i0, i1) get_block(i0, i1, 1L, nc),
    get_cols = function(j0, j1) get_block(1L, nr, j0, j1),
    get_block = get_block,
    storage_type = "filematrix",
    chunk_size = chunk_size
  )

  class(provider) <- c("filematrix_row_provider", "row_block_provider")
  provider
}

.validate_provider_dim <- function(x, label) {
  if (length(x) != 1L || is.na(x) || !is.finite(x) || x < 1L) {
    stop(label, " must be a positive finite dimension.", call. = FALSE)
  }
  invisible(TRUE)
}

.as_positive_scalar_integer <- function(x, label) {
  if (length(x) != 1L || is.na(x) || !is.finite(x) || x < 1L || x != as.integer(x)) {
    stop(label, " must be a positive scalar integer.", call. = FALSE)
  }
  as.integer(x)
}

.validate_closed_interval <- function(i0, i1, upper, label) {
  if (i0 > i1) {
    stop(label, " start must be <= ", label, " end.", call. = FALSE)
  }
  if (i1 > upper) {
    stop(label, " interval is outside provider dimensions.", call. = FALSE)
  }
  invisible(TRUE)
}

.check_filematrix_numeric_storage <- function(X) {
  storage_type <- tryCatch(X$type, error = function(e) NA_character_)
  if (is.null(storage_type) || !length(storage_type)) {
    storage_type <- NA_character_
  }
  storage_type <- as.character(storage_type)

  if (length(storage_type) && !is.na(storage_type[[1L]])) {
    if (!storage_type[[1L]] %in% c("double", "integer", "numeric")) {
      stop(
        "Only numeric filematrix storage is supported; got type '",
        storage_type[[1L]], "'.",
        call. = FALSE
      )
    }
    return(invisible(TRUE))
  }

  probe <- X[1L, 1L]
  if (!is.numeric(probe)) {
    stop("Only numeric filematrix storage is supported.", call. = FALSE)
  }
  invisible(TRUE)
}

.as_dense_double_matrix <- function(x, nrow, ncol) {
  if (is.matrix(x) && identical(dim(x), c(nrow, ncol))) {
    out <- x
  } else {
    out <- matrix(as.vector(x), nrow = nrow, ncol = ncol)
  }

  if (!is.numeric(out)) {
    stop("Provider block is not numeric; unsupported filematrix storage.", call. = FALSE)
  }

  storage.mode(out) <- "double"
  out
}
