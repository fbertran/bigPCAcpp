#include <Rcpp.h>
#include <R_ext/Lapack.h>
#include <algorithm>
#include <cmath>
#include <limits>
#include <vector>

using namespace Rcpp;

namespace {

double median_inplace(std::vector<double>& values) {
    const std::size_t n = values.size();
    if (n == 0) {
        return NA_REAL;
    }
    const std::size_t mid = n / 2;
    std::nth_element(values.begin(), values.begin() + mid, values.end());
    const double middle = values[mid];
    if (n % 2 == 1) {
        return middle;
    }
    std::nth_element(values.begin(), values.begin() + mid - 1, values.end());
    return 0.5 * (middle + values[mid - 1]);
}

double mad_about_zero(const std::vector<double>& values) {
    std::vector<double> abs_vals(values.size());
    for (std::size_t i = 0; i < values.size(); ++i) {
        abs_vals[i] = std::fabs(values[i]);
    }
    return median_inplace(abs_vals);
}

struct DenseSvdResult {
    NumericMatrix u;
    NumericVector d;
    NumericMatrix v;
};

DenseSvdResult dense_svd(const NumericMatrix& x, int nu, int nv) {
    const int m = x.nrow();
    const int n = x.ncol();
    const int min_dim = std::min(m, n);
    const int econ_cols = min_dim;

    NumericMatrix u_mat(m, econ_cols);
    NumericMatrix vt_mat(econ_cols, n);
    NumericMatrix a_work = clone(x);
    NumericVector singular_values(min_dim);

    const char jobz = 'S';
    const int lda = std::max(1, m);
    const int ldu = std::max(1, m);
    const int ldvt = std::max(1, econ_cols);

    bool used_dgesdd = false;
    int info = 0;

    int lwork = -1;
    double work_query = 0.0;
    std::vector<int> iwork(std::max(1, 8 * min_dim));
    F77_CALL(dgesdd)(&jobz,
                     &m,
                     &n,
                     a_work.begin(),
                     &lda,
                     singular_values.begin(),
                     u_mat.begin(),
                     &ldu,
                     vt_mat.begin(),
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
        a_work = clone(x);
        F77_CALL(dgesdd)(&jobz,
                         &m,
                         &n,
                         a_work.begin(),
                         &lda,
                         singular_values.begin(),
                         u_mat.begin(),
                         &ldu,
                         vt_mat.begin(),
                         &ldvt,
                         work.data(),
                         &lwork,
                         iwork.data(),
                         &info FCONE);
        if (info == 0) {
            used_dgesdd = true;
        }
    }

    if (!used_dgesdd) {
        const char jobu = 'S';
        const char jobvt = 'S';
        lwork = -1;
        work_query = 0.0;
        info = 0;
        a_work = clone(x);
        F77_CALL(dgesvd)(&jobu,
                         &jobvt,
                         &m,
                         &n,
                         a_work.begin(),
                         &lda,
                         singular_values.begin(),
                         u_mat.begin(),
                         &ldu,
                         vt_mat.begin(),
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
        a_work = clone(x);
        F77_CALL(dgesvd)(&jobu,
                         &jobvt,
                         &m,
                         &n,
                         a_work.begin(),
                         &lda,
                         singular_values.begin(),
                         u_mat.begin(),
                         &ldu,
                         vt_mat.begin(),
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

    NumericMatrix u_out(m, nu);
    for (int i = 0; i < m; ++i) {
        for (int k = 0; k < nu; ++k) {
            u_out(i, k) = u_mat(i, k);
        }
    }

    NumericMatrix v_out(n, nv);
    for (int j = 0; j < n; ++j) {
        for (int k = 0; k < nv; ++k) {
            v_out(j, k) = vt_mat(k, j);
        }
    }

    return DenseSvdResult{u_out, singular_values, v_out};
}

} // namespace

// [[Rcpp::export(name = ".svd_robust_cpp", rng = false)]]
List svd_robust_cpp(const NumericMatrix& x,
                    int ncomp,
                    int max_iter,
                    double tol,
                    double huber_k) {
    const int n = x.nrow();
    const int p = x.ncol();
    if (n == 0 || p == 0) {
        stop("`x` must contain at least one row and one column");
    }
    if (ncomp <= 0) {
        ncomp = std::min(n, p);
    }
    if (ncomp > std::min(n, p)) {
        ncomp = std::min(n, p);
    }
    if (max_iter <= 0) {
        stop("`max_iter` must be positive");
    }
    if (!std::isfinite(tol) || tol <= 0) {
        stop("`tol` must be a positive finite number");
    }
    if (!std::isfinite(huber_k) || huber_k <= 0) {
        stop("`huber_k` must be a positive finite number");
    }

    const int nu = std::min(n, ncomp);
    const int nv = std::min(p, ncomp);

    NumericVector weights(n, 1.0);
    NumericVector sqrt_weights(n, 1.0);
    std::vector<double> prev_singular(ncomp, 0.0);

    NumericMatrix weighted(n, p);
    for (int i = 0; i < n; ++i) {
        const double w = sqrt_weights[i];
        for (int j = 0; j < p; ++j) {
            weighted(i, j) = x(i, j) * w;
        }
    }

    DenseSvdResult svd_res = dense_svd(weighted, nu, nv);
    NumericMatrix last_u = svd_res.u;
    NumericMatrix last_v = svd_res.v;
    NumericVector last_d = svd_res.d;
    std::vector<double> singular(ncomp);
    std::vector<double> row_residuals(n);
    int iterations = 1;

    for (int iter = 1; iter <= max_iter; ++iter) {
        if (iter > 1) {
            for (int i = 0; i < n; ++i) {
                const double w = sqrt_weights[i];
                for (int j = 0; j < p; ++j) {
                    weighted(i, j) = x(i, j) * w;
                }
            }
            svd_res = dense_svd(weighted, nu, nv);
            last_u = svd_res.u;
            last_v = svd_res.v;
            last_d = svd_res.d;
            iterations = iter;
        }

        for (int k = 0; k < ncomp; ++k) {
            singular[static_cast<std::size_t>(k)] = last_d[k];
        }

        for (int i = 0; i < n; ++i) {
            const double inv_weight = (sqrt_weights[i] > 0) ? 1.0 / sqrt_weights[i] : 0.0;
            double row_sq = 0.0;
            for (int j = 0; j < p; ++j) {
                double fitted_weighted = 0.0;
                for (int k = 0; k < nv; ++k) {
                    fitted_weighted += last_u(i, k) * singular[static_cast<std::size_t>(k)] * last_v(j, k);
                }
                const double fitted = fitted_weighted * inv_weight;
                const double resid = x(i, j) - fitted;
                row_sq += resid * resid;
            }
            row_residuals[static_cast<std::size_t>(i)] = std::sqrt(row_sq);
        }

        double scale = mad_about_zero(row_residuals);
        if (!std::isfinite(scale) || scale <= std::numeric_limits<double>::epsilon()) {
            double sum_sq = 0.0;
            for (double value : row_residuals) {
                sum_sq += value * value;
            }
            if (n > 0) {
                scale = std::sqrt(sum_sq / static_cast<double>(n));
            }
        }
        if (!std::isfinite(scale) || scale <= std::numeric_limits<double>::epsilon()) {
            iterations = iter;
            break;
        }

        NumericVector new_weights(n);
        for (int i = 0; i < n; ++i) {
            const double scaled = row_residuals[static_cast<std::size_t>(i)] / (scale * huber_k);
            double w;
            if (!std::isfinite(scaled) || scaled <= 1.0) {
                w = 1.0;
            } else {
                w = 1.0 / scaled;
            }
            if (!std::isfinite(w)) {
                w = 1.0;
            }
            if (w < 1e-6) {
                w = 1e-6;
            } else if (w > 1.0) {
                w = 1.0;
            }
            new_weights[i] = w;
        }

        double weight_change = 0.0;
        for (int i = 0; i < n; ++i) {
            weight_change = std::max(weight_change, std::fabs(new_weights[i] - weights[i]));
        }
        if (!std::isfinite(weight_change)) {
            weight_change = 0.0;
        }

        double rel_change = std::numeric_limits<double>::infinity();
        if (iter > 1) {
            double max_prev = 0.0;
            for (double val : prev_singular) {
                if (val > max_prev) {
                    max_prev = val;
                }
            }
            if (max_prev < 1.0) {
                max_prev = 1.0;
            }
            double max_diff = 0.0;
            for (int k = 0; k < ncomp; ++k) {
                const double diff = std::fabs(singular[static_cast<std::size_t>(k)] - prev_singular[static_cast<std::size_t>(k)]);
                if (diff > max_diff) {
                    max_diff = diff;
                }
            }
            rel_change = max_diff / max_prev;
        }

        prev_singular = singular;
        weights = new_weights;
        for (int i = 0; i < n; ++i) {
            sqrt_weights[i] = std::sqrt(weights[i]);
        }

        if (weight_change < tol && rel_change < tol) {
            iterations = iter;
            break;
        }
        iterations = iter;
    }

    NumericMatrix u_out(n, nu);
    for (int i = 0; i < n; ++i) {
        for (int k = 0; k < nu; ++k) {
            u_out(i, k) = last_u(i, k);
        }
    }

    NumericMatrix v_out(p, nv);
    for (int j = 0; j < p; ++j) {
        for (int k = 0; k < nv; ++k) {
            v_out(j, k) = last_v(j, k);
        }
    }

    NumericVector d_out(ncomp);
    for (int k = 0; k < ncomp; ++k) {
        d_out[k] = last_d[k];
    }

    return List::create(
        Named("u") = u_out,
        Named("d") = d_out,
        Named("v") = v_out,
        Named("weights") = weights,
        Named("iterations") = iterations
    );
}
