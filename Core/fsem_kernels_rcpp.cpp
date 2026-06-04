// fsem_kernels_rcpp.cpp
//
// Additional Rcpp kernels used by em.estimation.rcpp() to push the entire
// per-(EM-iteration, factor/regression) CV + inner-fit pattern into C++.
//
// Replaces the R-level nested loops
//      for (l in 1:length(delt)) for (it in 1:no.it) for (m in 1:mont) { ... }
// with a single C++ call per factor/regression, eliminating thousands of R
// function-call traversals per EM iteration.
//
// Plain Rcpp + LAPACK only (no Armadillo). Matrices are assumed to be in R's
// column-major layout.

#include <Rcpp.h>
#include <R_ext/Lapack.h>
#include <R_ext/BLAS.h>
#include <vector>
#include <algorithm>
using namespace Rcpp;

// ---- LAPACK helpers ----------------------------------------------------

// In-place general solve: A X = B, A is n x n (overwritten by LU), B is n x nrhs
// (overwritten by solution). Returns LAPACK info code (0 = success).
static int solve_inplace_general(NumericMatrix& A, NumericMatrix& B) {
  int n = A.nrow();
  int nrhs = B.ncol();
  int info = 0;
  std::vector<int> ipiv(n);
  F77_CALL(dgesv)(&n, &nrhs, A.begin(), &n, ipiv.data(),
                  B.begin(), &n, &info);
  return info;
}

// Returns inv(M) (n x n), via dgesv on the identity. Does not modify input.
static NumericMatrix mat_inverse(const NumericMatrix& M) {
  int n = M.nrow();
  NumericMatrix A(n, n);
  std::copy(M.begin(), M.end(), A.begin());
  NumericMatrix B(n, n);
  for (int i = 0; i < n; ++i) B(i, i) = 1.0;
  int info = solve_inplace_general(A, B);
  if (info != 0) Rcpp::stop("fsem_kernels_rcpp: matrix inversion failed (info=%d)", info);
  return B;
}

// ---- core block ops ----------------------------------------------------

// Accumulate (1/M) * sum_k t(F_k) %*% S %*% F_k into out_A (p x p)
// and       (1/M) * sum_k t(F_k) %*% S %*% r_k into out_b (p).
// F is (n*nb) x p (column-major), r is n*nb.
static void accumulate_xtSx_xtSr(const NumericMatrix& F,
                                 const NumericMatrix& Si,
                                 const double* r,
                                 int nb,
                                 double invM,
                                 NumericMatrix& out_A,
                                 double* out_b) {
  const int nr = F.nrow();
  const int p  = F.ncol();
  const int n  = nr / nb;
  std::vector<double> SF((size_t)nb * (size_t)p);
  std::vector<double> Sr((size_t)nb);
  for (int k = 0; k < n; ++k) {
    const int off = k * nb;
    // SF = Si %*% F_k   (nb x p)
    for (int c = 0; c < p; ++c) {
      for (int rr = 0; rr < nb; ++rr) {
        double s = 0.0;
        for (int t = 0; t < nb; ++t) s += Si(rr, t) * F(off + t, c);
        SF[(size_t)rr + (size_t)c * nb] = s;
      }
    }
    // Sr = Si %*% r_k
    for (int rr = 0; rr < nb; ++rr) {
      double s = 0.0;
      for (int t = 0; t < nb; ++t) s += Si(rr, t) * r[off + t];
      Sr[rr] = s;
    }
    // out_b += invM * t(F_k) %*% Sr
    // out_A += invM * t(F_k) %*% SF
    for (int a = 0; a < p; ++a) {
      double sa = 0.0;
      for (int rr = 0; rr < nb; ++rr) sa += F(off + rr, a) * Sr[rr];
      out_b[a] += invM * sa;
      for (int b = 0; b < p; ++b) {
        double s = 0.0;
        for (int rr = 0; rr < nb; ++rr) s += F(off + rr, a) * SF[(size_t)rr + (size_t)b * nb];
        out_A(a, b) += invM * s;
      }
    }
  }
}

// Accumulate (1/(N*M)) * sum_k (r_k - F_k lambda)(r_k - F_k lambda)^T
// into out_sigma (nb x nb). F is (n*nb) x p, lambda is length p.
static void accumulate_sigma(const NumericMatrix& F,
                             const double* r,
                             const double* lambda,
                             int nb,
                             double invNM,
                             NumericMatrix& out_sigma) {
  const int nr = F.nrow();
  const int p  = F.ncol();
  const int n  = nr / nb;
  std::vector<double> ss((size_t)nb);
  for (int k = 0; k < n; ++k) {
    const int off = k * nb;
    for (int rr = 0; rr < nb; ++rr) {
      double s = r[off + rr];
      for (int c = 0; c < p; ++c) s -= F(off + rr, c) * lambda[c];
      ss[rr] = s;
    }
    for (int a = 0; a < nb; ++a) {
      const double sa = ss[a];
      for (int b = 0; b < nb; ++b)
        out_sigma(a, b) += invNM * sa * ss[b];
    }
  }
}

// ---- inner fit: replaces the for (it in 1:no.it) loop -------------------

// Inputs:
//   Fs    : List of length M of NumericMatrix, each (n_use * nb) x p
//   rs    : List of length M of NumericVector, each length n_use * nb
//           (these are z - f2 precomputed on R side)
//   sig_init : nb x nb starting residual covariance
//   pen      : p x p penalty matrix (the constant block, multiplied by delta)
//   delta    : scalar penalty weight
//   no_it    : number of inner iterations (default 10 in original)
//   nb       : basis size
//
// Output: List(lambda = NumericVector(p), sig = NumericMatrix(nb, nb))
//
// Exact equivalent of:
//   for (it in 1:no_it) {
//     Si <- solve(sig)
//     lam.sum  <- (1/M) sum_m xtSx(Fs[[m]], Si, nb) + delta*pen
//     lam1.sum <- (1/M) sum_m xtSz(Fs[[m]], Si, rs[[m]], nb)
//     lambda <- solve(lam.sum, lam1.sum)
//     sig <- (1/(n_use*M)) sum_m sum_k (rs_m,k - F_m,k lambda)(...)^T
//   }
// [[Rcpp::export]]
List fsem_fit_inner(List Fs, List rs, NumericMatrix sig_init,
                    NumericMatrix pen, double delta, int no_it, int nb) {
  const int M = Fs.size();
  if (M == 0) Rcpp::stop("fsem_fit_inner: Fs is empty");

  // Pre-extract pointers to avoid repeated as<>() in inner loops.
  std::vector<NumericMatrix> F_vec; F_vec.reserve(M);
  std::vector<NumericVector> r_vec; r_vec.reserve(M);
  for (int m = 0; m < M; ++m) {
    F_vec.push_back(as<NumericMatrix>(Fs[m]));
    r_vec.push_back(as<NumericVector>(rs[m]));
  }

  const int p        = F_vec[0].ncol();
  const int n_use_nb = F_vec[0].nrow();
  const int n_use    = n_use_nb / nb;
  const double invM  = 1.0 / (double)M;
  const double invNM = 1.0 / ((double)n_use * (double)M);

  NumericMatrix sig = clone(sig_init);
  NumericVector lambda(p);

  for (int it = 0; it < no_it; ++it) {
    NumericMatrix Si = mat_inverse(sig);

    // Build lam_sum (p x p) and lam1_sum (p)
    NumericMatrix lam_sum(p, p);
    std::vector<double> lam1_sum((size_t)p, 0.0);
    for (int m = 0; m < M; ++m) {
      accumulate_xtSx_xtSr(F_vec[m], Si, &(r_vec[m])[0], nb, invM, lam_sum, lam1_sum.data());
    }
    // Add delta * pen
    for (int a = 0; a < p; ++a)
      for (int b = 0; b < p; ++b)
        lam_sum(a, b) += delta * pen(a, b);

    // Solve lam_sum %*% lambda = lam1_sum
    NumericMatrix A = clone(lam_sum);
    NumericMatrix B(p, 1);
    for (int i = 0; i < p; ++i) B(i, 0) = lam1_sum[i];
    int info = solve_inplace_general(A, B);
    if (info != 0) Rcpp::stop("fsem_fit_inner: solve(lam_sum) failed (info=%d)", info);
    for (int i = 0; i < p; ++i) lambda[i] = B(i, 0);

    // sigma1 (nb x nb)
    NumericMatrix sigma1(nb, nb);
    for (int m = 0; m < M; ++m) {
      accumulate_sigma(F_vec[m], &(r_vec[m])[0], &lambda[0], nb, invNM, sigma1);
    }
    sig = sigma1;
  }

  return List::create(_["lambda"] = lambda, _["sig"] = sig);
}

// ---- full CV: replaces the for (l in 1:length(delt)) ... for (it ...) ---

// Runs the full CV path:
//   sig = sig_init                             (NOT reset between deltas - matches original)
//   for l in 1:D:
//     fit (no_it inner iterations on TRAIN) -> updates sig, computes lambda_l
//     Si_test = solve(sig)
//     per[l]  = (1/M) sum_m { -t(lambda_l) %*% b_m + t(lambda_l) %*% (A_m + delta*pen) %*% lambda_l }
// At the end: returns per, lambda from last delta, sig from last delta, and
// A_test/b_test caches from the last delta (used by R-side bound selection).
//
// Inputs:
//   Fs_train, rs_train : lists of length M (training fold)
//   Fs_test , rs_test  : lists of length M (test fold)
//   sig_init           : nb x nb
//   pen                : p x p
//   deltas             : numeric vector length D
//   no_it              : int
//   nb                 : int
//
// Output: List(per, lambda, sig, A_test, b_test).
// [[Rcpp::export]]
List fsem_cv_fit(List Fs_train, List rs_train,
                 List Fs_test,  List rs_test,
                 NumericMatrix sig_init, NumericMatrix pen,
                 NumericVector deltas, int no_it, int nb) {
  const int D = deltas.size();
  const int M = Fs_test.size();
  if (D == 0) Rcpp::stop("fsem_cv_fit: deltas is empty");
  if (M == 0) Rcpp::stop("fsem_cv_fit: Fs_test is empty");

  // Pre-extract test data pointers
  std::vector<NumericMatrix> Ft_vec; Ft_vec.reserve(M);
  std::vector<NumericVector> rt_vec; rt_vec.reserve(M);
  for (int m = 0; m < M; ++m) {
    Ft_vec.push_back(as<NumericMatrix>(Fs_test[m]));
    rt_vec.push_back(as<NumericVector>(rs_test[m]));
  }

  NumericMatrix sig = clone(sig_init);
  NumericVector per(D);
  NumericVector lambda_last;
  NumericMatrix sig_last;
  List A_test_last(M);
  List b_test_last(M);

  for (int l = 0; l < D; ++l) {
    const double delta = deltas[l];
    // Inner fit (mutates sig)
    List fit = fsem_fit_inner(Fs_train, rs_train, sig, pen, delta, no_it, nb);
    NumericVector lambda = fit["lambda"];
    sig = as<NumericMatrix>(fit["sig"]);
    const int p = lambda.size();

    // Compute test likelihood and (for the last delta) cache A_test/b_test
    NumericMatrix Si_test = mat_inverse(sig);
    double likeli = 0.0;
    const double invM = 1.0 / (double)M;

    for (int m = 0; m < M; ++m) {
      const NumericMatrix& Fm = Ft_vec[m];
      const NumericVector& rm = rt_vec[m];
      const int nr = Fm.nrow();
      const int n_use_test = nr / nb;

      NumericMatrix A(p, p);
      NumericVector b(p);
      std::vector<double> SF((size_t)nb * (size_t)p);
      std::vector<double> Sr((size_t)nb);
      for (int k = 0; k < n_use_test; ++k) {
        const int off = k * nb;
        for (int c = 0; c < p; ++c) {
          for (int rr = 0; rr < nb; ++rr) {
            double s = 0.0;
            for (int t = 0; t < nb; ++t) s += Si_test(rr, t) * Fm(off + t, c);
            SF[(size_t)rr + (size_t)c * nb] = s;
          }
        }
        for (int rr = 0; rr < nb; ++rr) {
          double s = 0.0;
          for (int t = 0; t < nb; ++t) s += Si_test(rr, t) * rm[off + t];
          Sr[rr] = s;
        }
        for (int a = 0; a < p; ++a) {
          double sa = 0.0;
          for (int rr = 0; rr < nb; ++rr) sa += Fm(off + rr, a) * Sr[rr];
          b[a] += sa;
          for (int b_ = 0; b_ < p; ++b_) {
            double s = 0.0;
            for (int rr = 0; rr < nb; ++rr) s += Fm(off + rr, a) * SF[(size_t)rr + (size_t)b_ * nb];
            A(a, b_) += s;
          }
        }
      }

      // val = -t(lambda) %*% b + t(lambda) %*% (A + delta*pen) %*% lambda
      double val = 0.0;
      for (int i = 0; i < p; ++i) val -= lambda[i] * b[i];
      for (int i = 0; i < p; ++i) {
        double s = 0.0;
        for (int j = 0; j < p; ++j) s += (A(i, j) + delta * pen(i, j)) * lambda[j];
        val += lambda[i] * s;
      }
      likeli += invM * val;

      if (l == D - 1) {
        A_test_last[m] = A;
        b_test_last[m] = b;
      }
    }
    per[l] = likeli;

    if (l == D - 1) {
      lambda_last = lambda;
      sig_last = sig;
    }
  }

  return List::create(_["per"] = per,
                      _["lambda"] = lambda_last,
                      _["sig"] = sig_last,
                      _["A_test"] = A_test_last,
                      _["b_test"] = b_test_last);
}
