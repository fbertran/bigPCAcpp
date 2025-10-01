#include <Rcpp.h>
#include <R_ext/BLAS.h>
#include <R_ext/Lapack.h>
#include <bigmemory/BigMatrix.h>
#include <bigmemory/MatrixAccessor.hpp>
#include <algorithm>
#include <cmath>
#include <limits>
#include <string>
#include <vector>

using namespace Rcpp;

namespace {

template <typename T>
List pca_impl(BigMatrix* pMat,
              bool center,
              bool scale,
              int ncomp,
              std::size_t block_size) {
  if (pMat == nullptr) {
    stop("Null big.matrix pointer");
  }
  if (pMat->ncol() == 0 || pMat->nrow() == 0) {
    stop("big.matrix must have at least one row and one column");
  }
  if (pMat->nrow() < 2) {
    stop("PCA requires at least two observations");
  }
  if (scale && !center) {
    stop("Scaling requires centering the variables");
  }
  if (block_size == 0) {
    stop("block_size must be strictly positive");
  }

  const std::size_t n = static_cast<std::size_t>(pMat->nrow());
  const std::size_t p = static_cast<std::size_t>(pMat->ncol());

  MatrixAccessor<T> accessor(*pMat);

  // Use Welford's algorithm for numerically stable column statistics.
  std::vector<long double> mean_acc(p, 0.0L);
  std::vector<long double> m2_acc(p, 0.0L);
  for (std::size_t row = 0; row < n; ++row) {
    const long double count = static_cast<long double>(row + 1);
    for (std::size_t col = 0; col < p; ++col) {
      const long double value = static_cast<long double>(accessor[col][row]);
      const long double delta = value - mean_acc[col];
      mean_acc[col] += delta / count;
      const long double delta2 = value - mean_acc[col];
      m2_acc[col] += delta * delta2;
    }
  }

  NumericVector centerVec(center ? static_cast<R_xlen_t>(p) : 0);
  if (center) {
    for (std::size_t j = 0; j < p; ++j) {
      centerVec[static_cast<R_xlen_t>(j)] = static_cast<double>(mean_acc[j]);
    }
  }

  NumericVector columnSd(static_cast<R_xlen_t>(p));
  NumericVector scaleVec(scale ? static_cast<R_xlen_t>(p) : 0);
  for (std::size_t j = 0; j < p; ++j) {
    const long double variance = (n > 1) ? (m2_acc[j] / static_cast<long double>(n - 1)) : 0.0L;
    const double sd = std::sqrt(std::max<double>(static_cast<double>(variance), 0.0));
    columnSd[static_cast<R_xlen_t>(j)] = sd;
    if (scale) {
      scaleVec[static_cast<R_xlen_t>(j)] = (sd == 0.0) ? 1.0 : sd;
    }
  }

  NumericMatrix covariance(static_cast<R_xlen_t>(p), static_cast<R_xlen_t>(p));
  const std::size_t effective_block = std::max<std::size_t>(1, std::min(block_size, n));
  std::vector<double> block_buffer(effective_block * p);
  bool first_block = true;

  for (std::size_t block_start = 0; block_start < n; block_start += effective_block) {
    const std::size_t block_end = std::min(block_start + effective_block, n);
    const std::size_t block_rows = block_end - block_start;

    // Populate column-major block with centered/scaled data.
    for (std::size_t col = 0; col < p; ++col) {
      const double center_val = center ? static_cast<double>(mean_acc[col]) : 0.0;
      const double scale_val = (scale ? scaleVec[static_cast<R_xlen_t>(col)] : 1.0);
      const double inv_scale = (scale_val == 0.0) ? 1.0 : 1.0 / scale_val;
      for (std::size_t row = 0; row < block_rows; ++row) {
        const std::size_t global_row = block_start + row;
        double value = static_cast<double>(accessor[col][global_row]);
        if (center) {
          value -= center_val;
        }
        value *= inv_scale;
        block_buffer[col * block_rows + row] = value;
      }
    }

    const char uplo = 'U';
    const char trans = 'T';
    const int n_int = static_cast<int>(p);
    const int k_int = static_cast<int>(block_rows);
    const double alpha = 1.0;
    const double beta = first_block ? 0.0 : 1.0;
    F77_CALL(dsyrk)(&uplo,
                   &trans,
                   &n_int,
                   &k_int,
                   &alpha,
                   block_buffer.data(),
                   &k_int,
                   &beta,
                   covariance.begin(),
                   &n_int FCONE FCONE);
    first_block = false;
  }

  const double denom = static_cast<double>(n - 1);
  if (denom <= 0.0) {
    stop("Unable to compute covariance with fewer than two observations");
  }
  double* cov_ptr = covariance.begin();
  for (std::size_t col = 0; col < p; ++col) {
    for (std::size_t col2 = col; col2 < p; ++col2) {
      const std::size_t idx = col + p * col2;
      const double value = cov_ptr[idx] / denom;
      cov_ptr[idx] = value;
      if (col2 != col) {
        cov_ptr[col2 + p * col] = value;
      }
    }
  }

  NumericMatrix eigenvectors = clone(covariance);
  std::vector<double> eigenvalues(p);

  char jobz = 'V';
  char uplo = 'U';
  int dim = static_cast<int>(p);
  int info = 0;
  int lwork = -1;
  double work_query = 0.0;
  F77_CALL(dsyev)(&jobz, &uplo, &dim, eigenvectors.begin(), &dim,
                  eigenvalues.data(), &work_query, &lwork, &info FCONE FCONE);
  if (info != 0) {
    stop("LAPACK dsyev workspace query failed");
  }
  lwork = static_cast<int>(work_query);
  if (lwork < 3 * dim) {
    lwork = 3 * dim;
  }
  std::vector<double> work(static_cast<std::size_t>(lwork));
  F77_CALL(dsyev)(&jobz, &uplo, &dim, eigenvectors.begin(), &dim,
                  eigenvalues.data(), work.data(), &lwork, &info FCONE FCONE);
  if (info != 0) {
    stop("LAPACK dsyev failed to converge");
  }

  int components_to_keep = ncomp;
  if (components_to_keep <= 0 || components_to_keep > dim) {
    components_to_keep = dim;
  }

  NumericMatrix rotation(static_cast<R_xlen_t>(p), components_to_keep);
  NumericVector eigenvalues_out(components_to_keep);
  NumericVector sdev(components_to_keep);
  double total_variance = 0.0;
  for (double v : eigenvalues) {
    total_variance += v;
  }
  NumericVector explained(components_to_keep);
  NumericVector cumulative(components_to_keep);

  for (int comp = 0; comp < components_to_keep; ++comp) {
    int source = dim - 1 - comp;
    double eigenvalue = eigenvalues[static_cast<std::size_t>(source)];
    eigenvalues_out[comp] = eigenvalue;
    double sd = std::sqrt(std::max(eigenvalue, 0.0));
    sdev[comp] = sd;
    for (int row = 0; row < dim; ++row) {
      rotation(row, comp) = eigenvectors(row, source);
    }
  }

  if (total_variance > 0.0) {
    double cumulative_sum = 0.0;
    for (int comp = 0; comp < components_to_keep; ++comp) {
      explained[comp] = eigenvalues_out[comp] / total_variance;
      cumulative_sum += explained[comp];
      cumulative[comp] = cumulative_sum;
    }
  }

  return List::create(
    Named("sdev") = sdev,
    Named("rotation") = rotation,
    Named("center") = center ? SEXP(centerVec) : R_NilValue,
    Named("scale") = scale ? SEXP(scaleVec) : R_NilValue,
    Named("column_sd") = columnSd,
    Named("eigenvalues") = eigenvalues_out,
    Named("explained_variance") = explained,
    Named("cumulative_variance") = cumulative,
    Named("covariance") = covariance,
    Named("nobs") = static_cast<double>(n)
  );
}

void orthonormalize_columns(std::vector<double>& mat,
                            std::size_t nrow,
                            int ncol) {
  for (int col = 0; col < ncol; ++col) {
    double* current = mat.data() + static_cast<std::size_t>(col) * nrow;
    for (int prev = 0; prev < col; ++prev) {
      double* previous = mat.data() + static_cast<std::size_t>(prev) * nrow;
      long double dot = 0.0L;
      for (std::size_t row = 0; row < nrow; ++row) {
        dot += static_cast<long double>(current[row]) * previous[row];
      }
      for (std::size_t row = 0; row < nrow; ++row) {
        current[row] -= static_cast<double>(dot * previous[row]);
      }
    }
    long double norm_sq = 0.0L;
    for (std::size_t row = 0; row < nrow; ++row) {
      norm_sq += static_cast<long double>(current[row]) * current[row];
    }
    if (norm_sq <= std::numeric_limits<long double>::min()) {
      for (std::size_t row = 0; row < nrow; ++row) {
        current[row] = 0.0;
      }
      if (nrow > 0) {
        current[static_cast<std::size_t>(col) % nrow] = 1.0;
      }
      continue;
    }
    const double inv_norm = 1.0 / std::sqrt(static_cast<double>(norm_sq));
    for (std::size_t row = 0; row < nrow; ++row) {
      current[row] *= inv_norm;
    }
  }
}

double compute_subspace_delta(const std::vector<double>& q_new,
                              const std::vector<double>& q_old,
                              std::size_t nrow,
                              int ncol,
                              std::vector<double>& workspace) {
  if (workspace.size() < static_cast<std::size_t>(ncol) * static_cast<std::size_t>(ncol)) {
    workspace.resize(static_cast<std::size_t>(ncol) * static_cast<std::size_t>(ncol));
  }
  const int rows = static_cast<int>(nrow);
  const int cols = ncol;
  const double alpha = 1.0;
  const double beta = 0.0;
  F77_CALL(dgemm)("T",
                  "N",
                  &cols,
                  &cols,
                  &rows,
                  &alpha,
                  const_cast<double*>(q_new.data()),
                  &rows,
                  const_cast<double*>(q_old.data()),
                  &rows,
                  &beta,
                  workspace.data(),
                  &cols FCONE FCONE);
  long double fro_sq = 0.0L;
  for (double value : workspace) {
    fro_sq += static_cast<long double>(value) * value;
  }
  const long double target = 2.0L * static_cast<long double>(cols) - 2.0L * fro_sq;
  const long double clamped = target < 0.0L ? 0.0L : target;
  return std::sqrt(static_cast<double>(clamped));
}

template <typename T>
List spca_power_impl(BigMatrix* pMat,
                     bool center,
                     bool scale,
                     int ncomp,
                     std::size_t block_size,
                     int max_iter,
                     double tol,
                     bool return_scores,
                     bool verbose) {
  if (pMat == nullptr) {
    stop("Null big.matrix pointer");
  }
  if (pMat->ncol() == 0 || pMat->nrow() == 0) {
    stop("big.matrix must have at least one row and one column");
  }
  if (pMat->nrow() < 2) {
    stop("PCA requires at least two observations");
  }
  if (scale && !center) {
    stop("Scaling requires centering the variables");
  }
  if (block_size == 0) {
    stop("block_size must be strictly positive");
  }
  if (max_iter <= 0) {
    stop("max_iter must be strictly positive");
  }
  if (tol < 0.0) {
    stop("tol must be non-negative");
  }

  const std::size_t n = static_cast<std::size_t>(pMat->nrow());
  const std::size_t p = static_cast<std::size_t>(pMat->ncol());

  int components_to_keep = ncomp;
  const int max_components = static_cast<int>(std::min<std::size_t>(n, p));
  if (components_to_keep <= 0 || components_to_keep > max_components) {
    components_to_keep = max_components;
  }
  if (components_to_keep <= 0) {
    stop("`ncomp` must select at least one component");
  }

  MatrixAccessor<T> accessor(*pMat);

  std::vector<long double> mean_acc(p, 0.0L);
  std::vector<long double> m2_acc(p, 0.0L);
  for (std::size_t row = 0; row < n; ++row) {
    const long double count = static_cast<long double>(row + 1);
    for (std::size_t col = 0; col < p; ++col) {
      const long double value = static_cast<long double>(accessor[col][row]);
      const long double delta = value - mean_acc[col];
      mean_acc[col] += delta / count;
      const long double delta2 = value - mean_acc[col];
      m2_acc[col] += delta * delta2;
    }
  }

  NumericVector centerVec(center ? static_cast<R_xlen_t>(p) : 0);
  if (center) {
    for (std::size_t j = 0; j < p; ++j) {
      centerVec[static_cast<R_xlen_t>(j)] = static_cast<double>(mean_acc[j]);
    }
  }

  NumericVector columnSd(static_cast<R_xlen_t>(p));
  NumericVector scaleVec(scale ? static_cast<R_xlen_t>(p) : 0);
  long double sumsq_acc = 0.0L;
  long double sumsq_scaled_acc = 0.0L;
  double total_variance = 0.0;
  for (std::size_t j = 0; j < p; ++j) {
    const long double variance = (n > 1)
      ? (m2_acc[j] / static_cast<long double>(n - 1))
      : 0.0L;
    const double sd = std::sqrt(std::max<double>(static_cast<double>(variance), 0.0));
    columnSd[static_cast<R_xlen_t>(j)] = sd;
    double scale_value = 1.0;
    if (scale) {
      scale_value = (sd == 0.0) ? 1.0 : sd;
      scaleVec[static_cast<R_xlen_t>(j)] = scale_value;
    }
    if (center) {
      total_variance += sd * sd;
    } else {
      const long double mean = mean_acc[j];
      const long double sumsq = m2_acc[j] + static_cast<long double>(n) * mean * mean;
      sumsq_acc += sumsq;
      if (scale) {
        const long double scale_sq = static_cast<long double>(scale_value) * scale_value;
        sumsq_scaled_acc += (scale_sq == 0.0L) ? 0.0L : (sumsq / scale_sq);
      }
    }
  }
  if (scale) {
    if (center) {
      total_variance = static_cast<double>(p);
    } else {
      total_variance = (n > 0)
        ? static_cast<double>(sumsq_scaled_acc / static_cast<long double>(n))
        : 0.0;
    }
  } else if (!center) {
    total_variance = (n > 0)
      ? static_cast<double>(sumsq_acc / static_cast<long double>(n))
      : 0.0;
  }

  const std::size_t effective_block = std::max<std::size_t>(1, std::min(block_size, n));
  std::vector<double> block_buffer(effective_block * p);
  std::vector<double> block_times_q(effective_block * static_cast<std::size_t>(components_to_keep));
  std::vector<double> operator_Q(p * static_cast<std::size_t>(components_to_keep));
  std::vector<double> delta_workspace(static_cast<std::size_t>(components_to_keep) * static_cast<std::size_t>(components_to_keep));

  NumericVector rand = Rcpp::rnorm(static_cast<R_xlen_t>(p) * components_to_keep);
  std::vector<double> Q_data(operator_Q.size());
  std::copy(rand.begin(), rand.end(), Q_data.begin());
  orthonormalize_columns(Q_data, p, components_to_keep);

  std::vector<double> center_buffer(p, 0.0);
  if (center) {
    std::copy(centerVec.begin(), centerVec.end(), center_buffer.begin());
  }
  std::vector<double> scale_buffer(p, 1.0);
  if (scale) {
    std::copy(scaleVec.begin(), scaleVec.end(), scale_buffer.begin());
    for (double& value : scale_buffer) {
      if (value == 0.0) {
        value = 1.0;
      }
    }
  }

  const int p_int = static_cast<int>(p);
  const int k_int = components_to_keep;
  const int max_iter_int = max_iter;
  const double alpha = 1.0;

  auto fill_block = [&](std::size_t block_start,
                        std::size_t block_rows) {
    for (std::size_t col = 0; col < p; ++col) {
      const double center_val = center ? center_buffer[col] : 0.0;
      const double scale_val = scale ? scale_buffer[col] : 1.0;
      const double inv_scale = (scale_val == 0.0) ? 1.0 : 1.0 / scale_val;
      for (std::size_t row = 0; row < block_rows; ++row) {
        const std::size_t global_row = block_start + row;
        double value = static_cast<double>(accessor[col][global_row]);
        if (center) {
          value -= center_val;
        }
        value *= inv_scale;
        block_buffer[col * block_rows + row] = value;
      }
    }
  };

  auto apply_operator = [&](const std::vector<double>& Q_mat,
                            std::vector<double>& output) {
    std::fill(output.begin(), output.end(), 0.0);
    for (std::size_t block_start = 0; block_start < n; block_start += effective_block) {
      const std::size_t block_end = std::min(block_start + effective_block, n);
      const std::size_t block_rows = block_end - block_start;
      fill_block(block_start, block_rows);
      const int block_rows_int = static_cast<int>(block_rows);
      const double beta_zero = 0.0;
      F77_CALL(dgemm)("N",
                      "N",
                      &block_rows_int,
                      &k_int,
                      &p_int,
                      &alpha,
                      block_buffer.data(),
                      &block_rows_int,
                      const_cast<double*>(Q_mat.data()),
                      &p_int,
                      &beta_zero,
                      block_times_q.data(),
                      &block_rows_int FCONE FCONE);
      const double beta_one = 1.0;
      F77_CALL(dgemm)("T",
                      "N",
                      &p_int,
                      &k_int,
                      &block_rows_int,
                      &alpha,
                      block_buffer.data(),
                      &block_rows_int,
                      block_times_q.data(),
                      &block_rows_int,
                      &beta_one,
                      output.data(),
                      &p_int FCONE FCONE);
    }
  };

  bool converged = false;
  double final_delta = std::numeric_limits<double>::quiet_NaN();
  int iterations = 0;
  std::vector<double> Q_new(operator_Q.size());

  for (int iter = 0; iter < max_iter_int; ++iter) {
    iterations = iter + 1;
    apply_operator(Q_data, operator_Q);
    std::vector<double> op_copy = operator_Q;
    std::vector<double> singular(static_cast<std::size_t>(k_int));
    std::vector<double> vt_dummy(1);
    char jobu = 'S';
    char jobvt = 'N';
    int m = p_int;
    int ncols = k_int;
    int lda = p_int;
    int ldu = p_int;
    int ldvt = 1;
    int info = 0;
    int lwork = -1;
    double work_query = 0.0;
    F77_CALL(dgesvd)(&jobu,
                     &jobvt,
                     &m,
                     &ncols,
                     op_copy.data(),
                     &lda,
                     singular.data(),
                     Q_new.data(),
                     &ldu,
                     vt_dummy.data(),
                     &ldvt,
                     &work_query,
                     &lwork,
                     &info FCONE FCONE);
    if (info != 0) {
      stop("LAPACK dgesvd workspace query failed");
    }
    lwork = std::max(3 * std::min(m, ncols), static_cast<int>(work_query));
    std::vector<double> work(static_cast<std::size_t>(lwork));
    F77_CALL(dgesvd)(&jobu,
                     &jobvt,
                     &m,
                     &ncols,
                     operator_Q.data(),
                     &lda,
                     singular.data(),
                     Q_new.data(),
                     &ldu,
                     vt_dummy.data(),
                     &ldvt,
                     work.data(),
                     &lwork,
                     &info FCONE FCONE);
    if (info != 0) {
      stop("LAPACK dgesvd failed to converge");
    }

    final_delta = compute_subspace_delta(Q_new, Q_data, p, k_int, delta_workspace);
    if (verbose) {
      Rcout << "Iteration " << iterations << ", subspace delta = " << final_delta << std::endl;
    }
    Q_data.swap(Q_new);
    if (final_delta <= tol) {
      converged = true;
      break;
    }
  }

  apply_operator(Q_data, operator_Q);

  const double denom = center ? static_cast<double>(n - 1) : static_cast<double>(n);
  if (denom <= 0.0) {
    stop("Unable to compute covariance with fewer than two observations");
  }

  std::vector<double> small_matrix(static_cast<std::size_t>(k_int) * static_cast<std::size_t>(k_int));
  const double beta_zero = 0.0;
  F77_CALL(dgemm)("T",
                  "N",
                  &k_int,
                  &k_int,
                  &p_int,
                  &alpha,
                  Q_data.data(),
                  &p_int,
                  operator_Q.data(),
                  &p_int,
                  &beta_zero,
                  small_matrix.data(),
                  &k_int FCONE FCONE);
  for (double& value : small_matrix) {
    value /= denom;
  }
  for (int col = 0; col < k_int; ++col) {
    for (int row = col + 1; row < k_int; ++row) {
      const double sym = 0.5 * (small_matrix[row + col * k_int] + small_matrix[col + row * k_int]);
      small_matrix[row + col * k_int] = sym;
      small_matrix[col + row * k_int] = sym;
    }
  }

  NumericVector eigenvalues(static_cast<R_xlen_t>(k_int));
  char jobz = 'V';
  char uplo = 'U';
  int info = 0;
  int lwork = -1;
  double work_query = 0.0;
  F77_CALL(dsyev)(&jobz,
                  &uplo,
                  &k_int,
                  small_matrix.data(),
                  &k_int,
                  eigenvalues.begin(),
                  &work_query,
                  &lwork,
                  &info FCONE FCONE);
  if (info != 0) {
    stop("LAPACK dsyev workspace query failed");
  }
  lwork = std::max(3 * k_int, static_cast<int>(work_query));
  std::vector<double> work(static_cast<std::size_t>(lwork));
  F77_CALL(dsyev)(&jobz,
                  &uplo,
                  &k_int,
                  small_matrix.data(),
                  &k_int,
                  eigenvalues.begin(),
                  work.data(),
                  &lwork,
                  &info FCONE FCONE);
  if (info != 0) {
    stop("LAPACK dsyev failed to converge");
  }

  std::vector<double> eigenvectors_sorted(static_cast<std::size_t>(k_int) * static_cast<std::size_t>(k_int));
  NumericVector eigenvalues_out(static_cast<R_xlen_t>(k_int));
  NumericVector sdev(static_cast<R_xlen_t>(k_int));
  for (int comp = 0; comp < k_int; ++comp) {
    const int source = k_int - 1 - comp;
    const double value = eigenvalues[source];
    eigenvalues_out[comp] = value;
    sdev[comp] = std::sqrt(std::max(value, 0.0));
    for (int row = 0; row < k_int; ++row) {
      eigenvectors_sorted[row + comp * k_int] = small_matrix[row + source * k_int];
    }
  }

  std::vector<double> rotation_data(static_cast<std::size_t>(p_int) * static_cast<std::size_t>(k_int));
  F77_CALL(dgemm)("N",
                  "N",
                  &p_int,
                  &k_int,
                  &k_int,
                  &alpha,
                  Q_data.data(),
                  &p_int,
                  eigenvectors_sorted.data(),
                  &k_int,
                  &beta_zero,
                  rotation_data.data(),
                  &p_int FCONE FCONE);

  NumericMatrix rotation(static_cast<R_xlen_t>(p), static_cast<R_xlen_t>(k_int));
  std::copy(rotation_data.begin(), rotation_data.end(), rotation.begin());

  NumericVector explained(static_cast<R_xlen_t>(k_int));
  NumericVector cumulative(static_cast<R_xlen_t>(k_int));
  double cumulative_sum = 0.0;
  for (int comp = 0; comp < k_int; ++comp) {
    const double value = eigenvalues_out[comp];
    const double explained_value = (total_variance > 0.0) ? (value / total_variance) : 0.0;
    explained[comp] = explained_value;
    cumulative_sum += explained_value;
    cumulative[comp] = cumulative_sum;
  }

  NumericMatrix scores_matrix;
  if (return_scores) {
    scores_matrix = NumericMatrix(static_cast<R_xlen_t>(n), static_cast<R_xlen_t>(k_int));
    std::vector<double> block_scores(effective_block * static_cast<std::size_t>(k_int));
    std::size_t offset = 0;
    for (std::size_t block_start = 0; block_start < n; block_start += effective_block) {
      const std::size_t block_end = std::min(block_start + effective_block, n);
      const std::size_t block_rows = block_end - block_start;
      fill_block(block_start, block_rows);
      const int block_rows_int = static_cast<int>(block_rows);
      F77_CALL(dgemm)("N",
                      "N",
                      &block_rows_int,
                      &k_int,
                      &p_int,
                      &alpha,
                      block_buffer.data(),
                      &block_rows_int,
                      rotation.begin(),
                      &p_int,
                      &beta_zero,
                      block_scores.data(),
                      &block_rows_int FCONE FCONE);
      for (int comp = 0; comp < k_int; ++comp) {
        for (std::size_t row = 0; row < block_rows; ++row) {
          scores_matrix(static_cast<R_xlen_t>(offset + row), comp) =
            block_scores[static_cast<std::size_t>(comp) * block_rows + row];
        }
      }
      offset += block_rows;
    }
  }

  List result = List::create(
    Named("sdev") = sdev,
    Named("rotation") = rotation,
    Named("center") = center ? SEXP(centerVec) : R_NilValue,
    Named("scale") = scale ? SEXP(scaleVec) : R_NilValue,
    Named("scores") = return_scores ? SEXP(scores_matrix) : R_NilValue,
    Named("column_sd") = columnSd,
    Named("eigenvalues") = eigenvalues_out,
    Named("explained_variance") = explained,
    Named("cumulative_variance") = cumulative,
    Named("covariance") = R_NilValue,
    Named("nobs") = static_cast<double>(n)
  );
  result.attr("iterations") = iterations;
  result.attr("tolerance") = tol;
  result.attr("converged") = converged;
  result.attr("delta") = final_delta;
  return result;
}

} // anonymous namespace

// [[Rcpp::export(name = ".pca_spca_bigmatrix")]]
SEXP pca_spca_bigmatrix(SEXP xpMat,
                        bool center = true,
                        bool scale = false,
                        int ncomp = -1,
                        std::size_t block_size = 2048,
                        int max_iter = 50,
                        double tol = 1e-4,
                        bool return_scores = false,
                        bool verbose = false) {
  Rcpp::XPtr<BigMatrix> pMat(xpMat);
  const int matrix_type = pMat->matrix_type();
  if (matrix_type == BigMatrix::FLOAT || matrix_type == static_cast<int>(sizeof(float))) {
    return spca_power_impl<float>(pMat, center, scale, ncomp, block_size, max_iter, tol, return_scores, verbose);
  }
  if (matrix_type == BigMatrix::DOUBLE || matrix_type == static_cast<int>(sizeof(double))) {
    return spca_power_impl<double>(pMat, center, scale, ncomp, block_size, max_iter, tol, return_scores, verbose);
  }
  stop("Unsupported big.matrix type. Use a double- or float-based big.matrix.");
}

// [[Rcpp::export(name = ".pca_bigmatrix")]]
SEXP pca_bigmatrix(SEXP xpMat,
                   bool center = true,
                   bool scale = false,
                   int ncomp = -1,
                   std::size_t block_size = 1024) {
  Rcpp::XPtr<BigMatrix> pMat(xpMat);
  const int matrix_type = pMat->matrix_type();
  if (matrix_type == BigMatrix::FLOAT || matrix_type == static_cast<int>(sizeof(float))) {
    return pca_impl<float>(pMat, center, scale, ncomp, block_size);
  }
  if (matrix_type == BigMatrix::DOUBLE || matrix_type == static_cast<int>(sizeof(double))) {
    return pca_impl<double>(pMat, center, scale, ncomp, block_size);
  }
  stop("Unsupported big.matrix type. Use a double- or float-based big.matrix.");
}

namespace {

struct SingularVectorRequest {
  int nu;
  int nv;
};

SingularVectorRequest normalize_singular_vector_request(int nu,
                                                        int nv,
                                                        int min_dim,
                                                        int m,
                                                        int ncols) {
  if (nu < 0) {
    nu = min_dim;
  }
  if (nv < 0) {
    nv = min_dim;
  }
  if (nu > min_dim) {
    stop("`nu` cannot exceed min(nrow, ncol)");
  }
  if (nv > min_dim) {
    stop("`nv` cannot exceed min(nrow, ncol)");
  }
  if (nu > m) {
    stop("`nu` cannot exceed the number of rows");
  }
  if (nv > ncols) {
    stop("`nv` cannot exceed the number of columns");
  }
  return {nu, nv};
}

template <typename T>
List svd_impl(BigMatrix* pMat,
              int nu,
              int nv,
              std::size_t block_size,
              bool prefer_dgesdd) {
  if (pMat == nullptr) {
    stop("Null big.matrix pointer");
  }
  const std::size_t n = static_cast<std::size_t>(pMat->nrow());
  const std::size_t p = static_cast<std::size_t>(pMat->ncol());
  if (n == 0 || p == 0) {
    stop("big.matrix must have at least one row and one column");
  }
  if (block_size == 0) {
    stop("block_size must be strictly positive");
  }

  int m = static_cast<int>(n);
  int ncols = static_cast<int>(p);
  int min_dim = std::min(m, ncols);
  if (min_dim <= 0) {
    stop("Unable to compute SVD with zero-sized dimension");
  }

  const SingularVectorRequest request = normalize_singular_vector_request(
    nu,
    nv,
    min_dim,
    m,
    ncols
  );
  nu = request.nu;
  nv = request.nv;

  MatrixAccessor<T> accessor(*pMat);
  NumericMatrix dense(m, ncols);

  const std::size_t effective_block = std::max<std::size_t>(1, std::min(block_size, n));
  std::vector<double> block_buffer(effective_block * p);
  double* dense_ptr = dense.begin();

  for (std::size_t block_start = 0; block_start < n; block_start += effective_block) {
    const std::size_t block_end = std::min(block_start + effective_block, n);
    const std::size_t block_rows = block_end - block_start;

    for (std::size_t col = 0; col < p; ++col) {
      for (std::size_t row = 0; row < block_rows; ++row) {
        const std::size_t global_row = block_start + row;
        block_buffer[col * block_rows + row] =
          static_cast<double>(accessor[col][global_row]);
      }
    }

    for (std::size_t col = 0; col < p; ++col) {
      double* dest = dense_ptr + static_cast<std::size_t>(col) * n + block_start;
      std::copy(block_buffer.data() + col * block_rows,
                block_buffer.data() + col * block_rows + block_rows,
                dest);
    }
  }

  const bool compute_u = nu > 0;
  const bool compute_v = nv > 0;
  char jobz = (compute_u || compute_v) ? 'S' : 'N';

  NumericMatrix u_mat;
  NumericMatrix vt_mat;
  std::vector<double> dummy_u;
  std::vector<double> dummy_vt;
  double* u_ptr = nullptr;
  double* vt_ptr = nullptr;
  int ldu = std::max(1, m);
  int ldvt = std::max(1, ncols);

  if (jobz == 'N') {
    dummy_u.resize(1);
    dummy_vt.resize(1);
    u_ptr = dummy_u.data();
    vt_ptr = dummy_vt.data();
    ldu = 1;
    ldvt = 1;
  } else { // 'S'
    const int econ_cols = min_dim;
    u_mat = NumericMatrix(m, econ_cols);
    vt_mat = NumericMatrix(econ_cols, ncols);
    u_ptr = u_mat.begin();
    vt_ptr = vt_mat.begin();
    ldu = std::max(1, m);
    ldvt = std::max(1, econ_cols);
  }

  NumericMatrix a_work = clone(dense);
  NumericVector singular_values(min_dim);

  bool used_dgesdd = false;
  int info = 0;

  if (prefer_dgesdd) {
    int lwork = -1;
    double work_query = 0.0;
    std::vector<int> iwork(std::max<int>(1, 8 * min_dim));
    F77_CALL(dgesdd)(&jobz,
                     &m,
                     &ncols,
                     a_work.begin(),
                     &m,
                     singular_values.begin(),
                     u_ptr,
                     &ldu,
                     vt_ptr,
                     &ldvt,
                     &work_query,
                     &lwork,
                     iwork.data(),
                     &info FCONE);
    if (info == 0) {
      lwork = static_cast<int>(work_query);
      if (lwork < 1) {
        lwork = 1;
      }
      std::vector<double> work(static_cast<std::size_t>(lwork));
      a_work = clone(dense);
      F77_CALL(dgesdd)(&jobz,
                       &m,
                       &ncols,
                       a_work.begin(),
                       &m,
                       singular_values.begin(),
                       u_ptr,
                       &ldu,
                       vt_ptr,
                       &ldvt,
                       work.data(),
                       &lwork,
                       iwork.data(),
                       &info FCONE);
      if (info == 0) {
        used_dgesdd = true;
      }
    }
  }

  if (!used_dgesdd) {
    char jobu = (jobz == 'A') ? 'A' : (jobz == 'S' ? 'S' : 'N');
    char jobvt = jobu;
    int lwork = -1;
    double work_query = 0.0;
    info = 0;
    a_work = clone(dense);
    F77_CALL(dgesvd)(&jobu,
                     &jobvt,
                     &m,
                     &ncols,
                     a_work.begin(),
                     &m,
                     singular_values.begin(),
                     u_ptr,
                     &ldu,
                     vt_ptr,
                     &ldvt,
                     &work_query,
                     &lwork,
                     &info FCONE FCONE);
    if (info != 0) {
      stop("LAPACK dgesvd workspace query failed");
    }
    lwork = static_cast<int>(work_query);
    if (lwork < std::max(1, 5 * min_dim)) {
      lwork = std::max(1, 5 * min_dim);
    }
    std::vector<double> work(static_cast<std::size_t>(lwork));
    a_work = clone(dense);
    F77_CALL(dgesvd)(&jobu,
                     &jobvt,
                     &m,
                     &ncols,
                     a_work.begin(),
                     &m,
                     singular_values.begin(),
                     u_ptr,
                     &ldu,
                     vt_ptr,
                     &ldvt,
                     work.data(),
                     &lwork,
                     &info FCONE FCONE);
    if (info < 0) {
      stop("Invalid argument supplied to LAPACK dgesvd");
    }
    if (info > 0) {
      stop("LAPACK dgesvd failed to converge");
    }
  } else {
    if (info < 0) {
      stop("Invalid argument supplied to LAPACK dgesdd");
    }
    if (info > 0) {
      stop("LAPACK dgesdd failed to converge");
    }
  }

  NumericVector d = clone(singular_values);

  NumericMatrix u_out(m, nu);
  if (nu > 0) {
    for (int col = 0; col < nu; ++col) {
      for (int row = 0; row < m; ++row) {
        u_out(row, col) = u_mat(row, col);
      }
    }
  }

  NumericMatrix v_out(ncols, nv);
  if (nv > 0) {
    for (int col = 0; col < nv; ++col) {
      for (int row = 0; row < ncols; ++row) {
        v_out(row, col) = vt_mat(col, row);
      }
    }
  }

  return List::create(
    Named("u") = u_out,
    Named("d") = d,
    Named("v") = v_out
  );
}

} // anonymous namespace

// [[Rcpp::export(name = ".svd_bigmatrix")]]
SEXP svd_bigmatrix(SEXP xpMat,
                   int nu = -1,
                   int nv = -1,
                   std::size_t block_size = 1024,
                   std::string method = "dgesdd") {
  Rcpp::XPtr<BigMatrix> pMat(xpMat);
  bool prefer_dgesdd;
  if (method == "dgesdd") {
    prefer_dgesdd = true;
  } else if (method == "dgesvd") {
    prefer_dgesdd = false;
  } else {
    stop("`method` must be either 'dgesdd' or 'dgesvd'");
  }
  const int matrix_type = pMat->matrix_type();
  if (matrix_type == BigMatrix::FLOAT || matrix_type == static_cast<int>(sizeof(float))) {
    return svd_impl<float>(pMat, nu, nv, block_size, prefer_dgesdd);
  }
  if (matrix_type == BigMatrix::DOUBLE || matrix_type == static_cast<int>(sizeof(double))) {
    return svd_impl<double>(pMat, nu, nv, block_size, prefer_dgesdd);
  }
  stop("Unsupported big.matrix type. Use a double- or float-based big.matrix.");
}

namespace {

template <typename T>
NumericMatrix scores_impl(BigMatrix* pMat,
                          const NumericMatrix& rotation,
                          const NumericVector& center,
                          const NumericVector& scale_vec,
                          int ncomp,
                          std::size_t block_size) {
  if (pMat == nullptr) {
    stop("Null big.matrix pointer");
  }
  if (block_size == 0) {
    stop("block_size must be strictly positive");
  }

  const std::size_t n = static_cast<std::size_t>(pMat->nrow());
  const std::size_t p = static_cast<std::size_t>(pMat->ncol());
  const int cols = rotation.ncol();
  if (rotation.nrow() != static_cast<int>(p)) {
    stop("Rotation matrix and big.matrix dimensions are incompatible");
  }

  int components_to_keep = ncomp;
  if (components_to_keep <= 0 || components_to_keep > cols) {
    components_to_keep = cols;
  }

  bool use_center = center.size() == static_cast<int>(p);
  bool use_scale = scale_vec.size() == static_cast<int>(p);
  if (use_scale && !use_center) {
    stop("Scaling requires providing column centers");
  }

  std::vector<double> center_buffer(p, 0.0);
  std::vector<double> scale_buffer(p, 1.0);
  if (use_center) {
    std::copy(center.begin(), center.end(), center_buffer.begin());
  }
  if (use_scale) {
    std::copy(scale_vec.begin(), scale_vec.end(), scale_buffer.begin());
    for (std::size_t j = 0; j < p; ++j) {
      if (scale_buffer[j] == 0.0) {
        scale_buffer[j] = 1.0;
      }
    }
  }

  MatrixAccessor<T> accessor(*pMat);
  NumericMatrix scores(static_cast<R_xlen_t>(n), components_to_keep);

  const std::size_t effective_block = std::max<std::size_t>(1, std::min(block_size, n));
  std::vector<double> block_buffer(effective_block * p);
  std::vector<double> block_scores(effective_block * static_cast<std::size_t>(components_to_keep));

  for (std::size_t block_start = 0; block_start < n; block_start += effective_block) {
    const std::size_t block_end = std::min(block_start + effective_block, n);
    const std::size_t block_rows = block_end - block_start;

    for (std::size_t col = 0; col < p; ++col) {
      const double center_val = use_center ? center_buffer[col] : 0.0;
      const double scale_val = use_scale ? scale_buffer[col] : 1.0;
      const double inv_scale = (scale_val == 0.0) ? 1.0 : 1.0 / scale_val;
      for (std::size_t row = 0; row < block_rows; ++row) {
        const std::size_t global_row = block_start + row;
        double value = static_cast<double>(accessor[col][global_row]);
        value -= center_val;
        value *= inv_scale;
        block_buffer[col * block_rows + row] = value;
      }
    }

    const char transA = 'N';
    const char transB = 'N';
    const int m = static_cast<int>(block_rows);
    const int n_comp = components_to_keep;
    const int k = static_cast<int>(p);
    const double alpha = 1.0;
    const double beta = 0.0;
    F77_CALL(dgemm)(&transA,
                    &transB,
                    &m,
                    &n_comp,
                    &k,
                    &alpha,
                    block_buffer.data(),
                    &m,
                    rotation.begin(),
                    &k,
                    &beta,
                    block_scores.data(),
                    &m FCONE FCONE);

    for (int comp = 0; comp < n_comp; ++comp) {
      for (std::size_t row = 0; row < block_rows; ++row) {
        scores(static_cast<int>(block_start + row), comp) =
          block_scores[static_cast<std::size_t>(comp) * block_rows + row];
      }
    }
  }

  return scores;
}

} // anonymous namespace

// [[Rcpp::export(name = ".pca_scores_bigmatrix")]]
SEXP pca_scores_bigmatrix(SEXP xpMat,
                          const NumericMatrix& rotation,
                          const NumericVector& center,
                          const NumericVector& scale,
                          int ncomp = -1,
                          std::size_t block_size = 1024) {
  Rcpp::XPtr<BigMatrix> pMat(xpMat);
  const int matrix_type = pMat->matrix_type();
  if (matrix_type == BigMatrix::FLOAT || matrix_type == static_cast<int>(sizeof(float))) {
    return scores_impl<float>(pMat, rotation, center, scale, ncomp, block_size);
  }
  if (matrix_type == BigMatrix::DOUBLE || matrix_type == static_cast<int>(sizeof(double))) {
    return scores_impl<double>(pMat, rotation, center, scale, ncomp, block_size);
  }
  stop("Unsupported big.matrix type. Use a double- or float-based big.matrix.");
}

namespace {

NumericMatrix compute_loadings(const NumericMatrix& rotation,
                               const NumericVector& sdev,
                               const NumericVector& scale_vec,
                               bool use_scale) {
  const int p = rotation.nrow();
  const int k = rotation.ncol();
  if (sdev.size() < k) {
    stop("Length of sdev must be at least the number of components");
  }
  if (use_scale && scale_vec.size() != p) {
    stop("Length of scale vector must match the number of variables");
  }

  NumericMatrix loadings(p, k);
  for (int comp = 0; comp < k; ++comp) {
    double sd = sdev[comp];
    for (int var = 0; var < p; ++var) {
      double value = rotation(var, comp) * sd;
      if (use_scale) {
        const double denom = scale_vec[var];
        if (std::abs(denom) > std::numeric_limits<double>::min()) {
          value /= denom;
        }
      }
      loadings(var, comp) = value;
    }
  }
  return loadings;
}

} // anonymous namespace

// [[Rcpp::export(name = ".pca_variable_loadings")]]
NumericMatrix pca_variable_loadings(const NumericMatrix& rotation,
                                    const NumericVector& sdev) {
  NumericVector empty_scale;
  return compute_loadings(rotation, sdev, empty_scale, false);
}

// [[Rcpp::export(name = ".pca_variable_correlations")]]
NumericMatrix pca_variable_correlations(const NumericMatrix& rotation,
                                        const NumericVector& sdev,
                                        const NumericVector& column_sd) {
  return compute_loadings(rotation, sdev, column_sd, true);
}

// [[Rcpp::export(name = ".pca_variable_contributions")]]
NumericMatrix pca_variable_contributions(const NumericMatrix& loadings) {
  const int p = loadings.nrow();
  const int k = loadings.ncol();
  NumericMatrix contributions(p, k);
  for (int comp = 0; comp < k; ++comp) {
    double denom = 0.0;
    for (int var = 0; var < p; ++var) {
      double weight = loadings(var, comp);
      denom += weight * weight;
    }
    if (denom <= std::numeric_limits<double>::min()) {
      continue;
    }
    for (int var = 0; var < p; ++var) {
      double weight = loadings(var, comp);
      contributions(var, comp) = (weight * weight) / denom;
    }
  }
  return contributions;
}

// [[Rcpp::export(name = ".pca_individual_contributions")]]
NumericMatrix pca_individual_contributions(const NumericMatrix& scores,
                                           const NumericVector& sdev,
                                           double total_weight = NA_REAL) {
  const int n = scores.nrow();
  const int k = scores.ncol();
  if (sdev.size() < k) {
    stop("Length of sdev must be at least the number of components");
  }
  double weight;
  if (R_IsNA(total_weight)) {
    if (n < 2) {
      stop("`scores` must contain at least two rows when `total_weight` is not provided");
    }
    weight = static_cast<double>(n - 1);
  } else {
    weight = total_weight;
  }
  if (!(weight > 0.0)) {
    stop("`total_weight` must be strictly positive");
  }

  NumericMatrix contributions(n, k);
  for (int comp = 0; comp < k; ++comp) {
    const double eigenvalue = sdev[comp] * sdev[comp];
    if (eigenvalue <= std::numeric_limits<double>::min()) {
      continue;
    }
    const double denom = weight * eigenvalue;
    for (int row = 0; row < n; ++row) {
      const double score = scores(row, comp);
      contributions(row, comp) = (score * score) / denom;
    }
  }
  return contributions;
}

// [[Rcpp::export(name = ".pca_individual_cos2")]]
NumericMatrix pca_individual_cos2(const NumericMatrix& scores) {
  const int n = scores.nrow();
  const int k = scores.ncol();
  NumericMatrix cos2(n, k);
  for (int row = 0; row < n; ++row) {
    double denom = 0.0;
    for (int comp = 0; comp < k; ++comp) {
      const double value = scores(row, comp);
      denom += value * value;
    }
    if (denom <= std::numeric_limits<double>::min()) {
      continue;
    }
    for (int comp = 0; comp < k; ++comp) {
      const double value = scores(row, comp);
      cos2(row, comp) = (value * value) / denom;
    }
  }
  return cos2;
}

// [[Rcpp::export(name = ".pca_variable_cos2")]]
NumericMatrix pca_variable_cos2(const NumericMatrix& correlations) {
  const int p = correlations.nrow();
  const int k = correlations.ncol();
  NumericMatrix cos2(p, k);
  for (int var = 0; var < p; ++var) {
    for (int comp = 0; comp < k; ++comp) {
      const double corr = correlations(var, comp);
      cos2(var, comp) = corr * corr;
    }
  }
  return cos2;
}

namespace {

inline void ensure_double_matrix(BigMatrix* ptr,
                                 const char* name,
                                 std::size_t expected_rows,
                                 std::size_t expected_cols) {
  if (ptr == nullptr) {
    stop("Null big.matrix pointer for %s", name);
  }
  if (ptr->matrix_type() != 8) {
    stop("%s big.matrix must store doubles", name);
  }
  if (ptr->nrow() != expected_rows || ptr->ncol() != expected_cols) {
    stop("%s big.matrix has incompatible dimensions", name);
  }
}

void copy_numeric_to_bigmatrix(const NumericMatrix& source,
                               BigMatrix* dest,
                               const char* name) {
  ensure_double_matrix(dest,
                       name,
                       static_cast<std::size_t>(source.nrow()),
                       static_cast<std::size_t>(source.ncol()));
  MatrixAccessor<double> dest_accessor(*dest);
  const int rows = source.nrow();
  const int cols = source.ncol();
  for (int col = 0; col < cols; ++col) {
    for (int row = 0; row < rows; ++row) {
      dest_accessor[col][row] = source(row, col);
    }
  }
}

template <typename T>
void scores_store_impl(BigMatrix* pMat,
                       BigMatrix* pDest,
                       const NumericMatrix& rotation,
                       const NumericVector& center,
                       const NumericVector& scale_vec,
                       int ncomp,
                       std::size_t block_size) {
  if (pMat == nullptr) {
    stop("Null big.matrix pointer");
  }
  if (pDest == nullptr) {
    stop("Null destination big.matrix pointer");
  }
  if (block_size == 0) {
    stop("block_size must be strictly positive");
  }

  const std::size_t n = static_cast<std::size_t>(pMat->nrow());
  const std::size_t p = static_cast<std::size_t>(pMat->ncol());
  const int cols = rotation.ncol();
  if (rotation.nrow() != static_cast<int>(p)) {
    stop("Rotation matrix and big.matrix dimensions are incompatible");
  }

  int components_to_keep = ncomp;
  if (components_to_keep <= 0 || components_to_keep > cols) {
    components_to_keep = cols;
  }

  ensure_double_matrix(pDest, "scores", n, static_cast<std::size_t>(components_to_keep));

  bool use_center = center.size() == static_cast<int>(p);
  bool use_scale = scale_vec.size() == static_cast<int>(p);
  if (use_scale && !use_center) {
    stop("Scaling requires providing column centers");
  }

  std::vector<double> center_buffer(p, 0.0);
  std::vector<double> scale_buffer(p, 1.0);
  if (use_center) {
    std::copy(center.begin(), center.end(), center_buffer.begin());
  }
  if (use_scale) {
    std::copy(scale_vec.begin(), scale_vec.end(), scale_buffer.begin());
    for (std::size_t j = 0; j < p; ++j) {
      if (scale_buffer[j] == 0.0) {
        scale_buffer[j] = 1.0;
      }
    }
  }

  MatrixAccessor<T> accessor(*pMat);
  MatrixAccessor<double> dest_accessor(*pDest);

  const std::size_t effective_block = std::max<std::size_t>(1, std::min(block_size, n));
  std::vector<double> block_buffer(effective_block * p);
  std::vector<double> block_scores(effective_block * static_cast<std::size_t>(components_to_keep));

  for (std::size_t block_start = 0; block_start < n; block_start += effective_block) {
    const std::size_t block_end = std::min(block_start + effective_block, n);
    const std::size_t block_rows = block_end - block_start;

    for (std::size_t col = 0; col < p; ++col) {
      const double center_val = use_center ? center_buffer[col] : 0.0;
      const double scale_val = use_scale ? scale_buffer[col] : 1.0;
      const double inv_scale = (scale_val == 0.0) ? 1.0 : 1.0 / scale_val;
      for (std::size_t row = 0; row < block_rows; ++row) {
        const std::size_t global_row = block_start + row;
        double value = static_cast<double>(accessor[col][global_row]);
        value -= center_val;
        value *= inv_scale;
        block_buffer[col * block_rows + row] = value;
      }
    }

    const char transA = 'N';
    const char transB = 'N';
    const int m = static_cast<int>(block_rows);
    const int n_comp = components_to_keep;
    const int k = static_cast<int>(p);
    const double alpha = 1.0;
    const double beta = 0.0;
    F77_CALL(dgemm)(&transA,
                    &transB,
                    &m,
                    &n_comp,
                    &k,
                    &alpha,
                    block_buffer.data(),
                    &m,
                    rotation.begin(),
                    &k,
                    &beta,
                    block_scores.data(),
                    &m FCONE FCONE);

    for (int comp = 0; comp < n_comp; ++comp) {
      for (std::size_t row = 0; row < block_rows; ++row) {
        dest_accessor[comp][block_start + row] =
          block_scores[static_cast<std::size_t>(comp) * block_rows + row];
      }
    }
  }
}

} // anonymous namespace

// [[Rcpp::export(name = ".pca_spca_stream_bigmatrix")]]
SEXP pca_spca_stream_bigmatrix(SEXP xpMat,
                               SEXP xpRotation,
                               bool center = true,
                               bool scale = false,
                               int ncomp = -1,
                               std::size_t block_size = 2048,
                               int max_iter = 50,
                               double tol = 1e-4,
                               bool return_scores = false,
                               bool verbose = false) {
  SEXP base_result = pca_spca_bigmatrix(xpMat, center, scale, ncomp, block_size, max_iter, tol, return_scores, verbose);
  if (Rf_isNull(xpRotation)) {
    return base_result;
  }
  Rcpp::List result(base_result);
  NumericMatrix rotation = result["rotation"];
  Rcpp::XPtr<BigMatrix> rotation_ptr(xpRotation);
  copy_numeric_to_bigmatrix(rotation, rotation_ptr, "rotation");
  result["rotation_stream_bigmatrix"] = xpRotation;
  return result;
}

// [[Rcpp::export(name = ".pca_stream_bigmatrix")]]
SEXP pca_stream_bigmatrix(SEXP xpMat,
                          SEXP xpRotation,
                          bool center = true,
                          bool scale = false,
                          int ncomp = -1,
                          std::size_t block_size = 1024) {
  SEXP base_result = pca_bigmatrix(xpMat, center, scale, ncomp, block_size);
  if (Rf_isNull(xpRotation)) {
    return base_result;
  }

  Rcpp::List result(base_result);
  NumericMatrix rotation = result["rotation"];
  Rcpp::XPtr<BigMatrix> rotation_ptr(xpRotation);
  copy_numeric_to_bigmatrix(rotation, rotation_ptr, "rotation");
  result["rotation_stream_bigmatrix"] = xpRotation;
  return result;
}

// [[Rcpp::export(name = ".pca_scores_stream_bigmatrix")]]
SEXP pca_scores_stream_bigmatrix(SEXP xpMat,
                                 SEXP xpDest,
                                 const NumericMatrix& rotation,
                                 const NumericVector& center,
                                 const NumericVector& scale,
                                 int ncomp = -1,
                                 std::size_t block_size = 1024) {
  Rcpp::XPtr<BigMatrix> pMat(xpMat);
  Rcpp::XPtr<BigMatrix> pDest(xpDest);
  const int matrix_type = pMat->matrix_type();
  if (matrix_type == BigMatrix::FLOAT || matrix_type == static_cast<int>(sizeof(float))) {
    scores_store_impl<float>(pMat, pDest, rotation, center, scale, ncomp, block_size);
  } else if (matrix_type == BigMatrix::DOUBLE || matrix_type == static_cast<int>(sizeof(double))) {
    scores_store_impl<double>(pMat, pDest, rotation, center, scale, ncomp, block_size);
  } else {
    stop("Unsupported big.matrix type. Use a double- or float-based big.matrix.");
  }
  return xpDest;
}

// [[Rcpp::export(name = ".pca_variable_loadings_stream_bigmatrix")]]
SEXP pca_variable_loadings_stream_bigmatrix(SEXP xpRotation,
                                            const NumericVector& sdev,
                                            SEXP xpDest) {
  Rcpp::XPtr<BigMatrix> rotation_ptr(xpRotation);
  Rcpp::XPtr<BigMatrix> dest_ptr(xpDest);
  ensure_double_matrix(rotation_ptr, "rotation", rotation_ptr->nrow(), rotation_ptr->ncol());
  const std::size_t p = static_cast<std::size_t>(rotation_ptr->nrow());
  const std::size_t k = static_cast<std::size_t>(rotation_ptr->ncol());
  if (sdev.size() < static_cast<int>(k)) {
    stop("Length of sdev must be at least the number of components");
  }
  ensure_double_matrix(dest_ptr, "loadings", p, k);

  MatrixAccessor<double> rotation_accessor(*rotation_ptr);
  MatrixAccessor<double> dest_accessor(*dest_ptr);
  for (std::size_t comp = 0; comp < k; ++comp) {
    double sd = sdev[static_cast<R_xlen_t>(comp)];
    for (std::size_t var = 0; var < p; ++var) {
      dest_accessor[comp][var] = rotation_accessor[comp][var] * sd;
    }
  }
  return xpDest;
}

// [[Rcpp::export(name = ".pca_variable_correlations_stream_bigmatrix")]]
SEXP pca_variable_correlations_stream_bigmatrix(SEXP xpRotation,
                                                const NumericVector& sdev,
                                                const NumericVector& column_sd,
                                                SEXP xpDest) {
  Rcpp::XPtr<BigMatrix> rotation_ptr(xpRotation);
  Rcpp::XPtr<BigMatrix> dest_ptr(xpDest);
  ensure_double_matrix(rotation_ptr, "rotation", rotation_ptr->nrow(), rotation_ptr->ncol());
  const std::size_t p = static_cast<std::size_t>(rotation_ptr->nrow());
  const std::size_t k = static_cast<std::size_t>(rotation_ptr->ncol());
  if (sdev.size() < static_cast<int>(k)) {
    stop("Length of sdev must be at least the number of components");
  }
  if (column_sd.size() != static_cast<int>(p)) {
    stop("Length of column_sd must match number of variables");
  }
  ensure_double_matrix(dest_ptr, "correlations", p, k);

  MatrixAccessor<double> rotation_accessor(*rotation_ptr);
  MatrixAccessor<double> dest_accessor(*dest_ptr);
  for (std::size_t comp = 0; comp < k; ++comp) {
    double sd = sdev[static_cast<R_xlen_t>(comp)];
    for (std::size_t var = 0; var < p; ++var) {
      double scale_val = column_sd[static_cast<R_xlen_t>(var)];
      if (std::abs(scale_val) <= std::numeric_limits<double>::min()) {
        dest_accessor[comp][var] = 0.0;
      } else {
        dest_accessor[comp][var] = (rotation_accessor[comp][var] * sd) / scale_val;
      }
    }
  }
  return xpDest;
}

// [[Rcpp::export(name = ".pca_variable_contributions_stream_bigmatrix")]]
SEXP pca_variable_contributions_stream_bigmatrix(SEXP xpLoadings,
                                                 SEXP xpDest) {
  Rcpp::XPtr<BigMatrix> loadings_ptr(xpLoadings);
  Rcpp::XPtr<BigMatrix> dest_ptr(xpDest);
  ensure_double_matrix(loadings_ptr, "loadings", loadings_ptr->nrow(), loadings_ptr->ncol());
  const std::size_t p = static_cast<std::size_t>(loadings_ptr->nrow());
  const std::size_t k = static_cast<std::size_t>(loadings_ptr->ncol());
  ensure_double_matrix(dest_ptr, "contributions", p, k);

  MatrixAccessor<double> load_accessor(*loadings_ptr);
  MatrixAccessor<double> dest_accessor(*dest_ptr);
  for (std::size_t comp = 0; comp < k; ++comp) {
    long double denom = 0.0L;
    for (std::size_t var = 0; var < p; ++var) {
      const long double weight = static_cast<long double>(load_accessor[comp][var]);
      denom += weight * weight;
    }
    if (denom <= std::numeric_limits<long double>::min()) {
      for (std::size_t var = 0; var < p; ++var) {
        dest_accessor[comp][var] = 0.0;
      }
      continue;
    }
    for (std::size_t var = 0; var < p; ++var) {
      const long double weight = static_cast<long double>(load_accessor[comp][var]);
      dest_accessor[comp][var] = static_cast<double>((weight * weight) / denom);
    }
  }
  return xpDest;
}
