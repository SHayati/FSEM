// fsem_kernels.cpp
// Optional Rcpp acceleration kernels for em.estimation.optimised().
//
// These kernels compute block quadratic / bilinear forms of the shape
//      t(F) %*% kron(I_n, S) %*% F      and      t(F) %*% kron(I_n, S) %*% z
// WITHOUT ever materialising the large (n*nb) x (n*nb) block-diagonal matrix
// kron(I_n, S). In the original code S = solve(sig) and the equivalent
// expression was solve(kronecker(diag(n), sig)) which is O((n*nb)^3); here it
// is reduced to n small (nb x nb) block operations.
//
// Plain Rcpp only (no Armadillo) so it compiles with Rcpp >= 1.0.

#include <Rcpp.h>
using namespace Rcpp;

// sum_k  t(F_k) %*% S %*% F_k     (F is (n*nb) x p, stacked in blocks of nb rows)
// [[Rcpp::export]]
NumericMatrix fsem_block_xtSx(const NumericMatrix& F, const NumericMatrix& S, int nb) {
  const int nr = F.nrow();
  const int p  = F.ncol();
  const int n  = nr / nb;
  NumericMatrix out(p, p);
  std::vector<double> SF((size_t)nb * p);
  for (int k = 0; k < n; ++k) {
    const int off = k * nb;
    // SF = S %*% F_k   (nb x p)
    for (int c = 0; c < p; ++c) {
      for (int r = 0; r < nb; ++r) {
        double s = 0.0;
        for (int t = 0; t < nb; ++t) s += S(r, t) * F(off + t, c);
        SF[(size_t)r + (size_t)c * nb] = s;
      }
    }
    // out += t(F_k) %*% SF
    for (int a = 0; a < p; ++a) {
      for (int b = 0; b < p; ++b) {
        double s = 0.0;
        for (int r = 0; r < nb; ++r) s += F(off + r, a) * SF[(size_t)r + (size_t)b * nb];
        out(a, b) += s;
      }
    }
  }
  return out;
}

// sum_k  t(F_k) %*% S %*% z_k      (z length n*nb)
// [[Rcpp::export]]
NumericVector fsem_block_xtSz(const NumericMatrix& F, const NumericMatrix& S,
                              const NumericVector& z, int nb) {
  const int nr = F.nrow();
  const int p  = F.ncol();
  const int n  = nr / nb;
  NumericVector out(p);
  std::vector<double> Sz(nb);
  for (int k = 0; k < n; ++k) {
    const int off = k * nb;
    for (int r = 0; r < nb; ++r) {
      double s = 0.0;
      for (int t = 0; t < nb; ++t) s += S(r, t) * z[off + t];
      Sz[r] = s;
    }
    for (int a = 0; a < p; ++a) {
      double s = 0.0;
      for (int r = 0; r < nb; ++r) s += F(off + r, a) * Sz[r];
      out[a] += s;
    }
  }
  return out;
}
