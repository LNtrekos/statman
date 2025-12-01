#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::plugins(cpp11)]]
// [[Rcpp::depends(Rcpp)]]

// ---------- Frank copula CDF ----------
inline double frank_copula(double u, double v, double theta) {

  if (std::abs(theta) < 1e-10)
    return u * v;  // independence

  double num = (std::exp(-theta * u) - 1.0) * (std::exp(-theta * v) - 1.0);
  double den = std::exp(-theta) - 1.0;

  return -(1.0 / theta) * std::log1p(num / den);

}

// [[Rcpp::export]]
NumericVector dmass_skellam2_fcop_cpp(
    NumericVector Z1,
    NumericVector Z2,
    double mu1, double mu2,
    double sigma1, double sigma2,
    double theta,
    Function pskellam2   // R Skellam CDF
) {
  int n = Z1.size();
  NumericVector out(n);

  for (int i = 0; i < n; i++) {

    double z1  = Z1[i];
    double z2  = Z2[i];

    // --- CDFs for Skellam2 ---
    double U1  = as<double>(pskellam2(z1,     mu1, sigma1));
    double U1m = as<double>(pskellam2(z1 - 1, mu1, sigma1));

    double V1  = as<double>(pskellam2(z2,     mu2, sigma2));
    double V1m = as<double>(pskellam2(z2 - 1, mu2, sigma2));

    // --- Bivariate mass via inclusion-exclusion ---
    double p =
      frank_copula(U1,  V1,  theta) -
      frank_copula(U1m, V1,  theta) -
      frank_copula(U1,  V1m, theta) +
      frank_copula(U1m, V1m, theta);

    if (p <= 0.0) p = 1e-16;

    out[i] = p;
  }

  return out;
}
