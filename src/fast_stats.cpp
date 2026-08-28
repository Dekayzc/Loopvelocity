// Fast primitives for per-TAD steady-state fits (Eq 13/14) and kNN metrics.
#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp;

//' Weighted least squares fit of the steady-state constraint
//'
//' Closed-form weighted linear regression y ~ a + b*x used for the per-TAD
//' steady-state fits (Eq 13/14). Weights are 0/1 selections of top/bottom
//' quantile cells.
//' @param x,y numeric vectors (I and S boundary contacts)
//' @param w nonnegative weights
//' @return list with \code{intercept} and \code{slope}
//' @examples
//' x <- 1:10; y <- 2 * x + 1
//' lv_wls_fit(x, y, rep(1, 10))
//' @export
// [[Rcpp::export]]
List lv_wls_fit(NumericVector x, NumericVector y, NumericVector w) {
  double sw = 0, swx = 0, swy = 0, swxx = 0, swxy = 0;
  for (int i = 0; i < x.size(); i++) {
    double wi = w[i];
    if (wi <= 0) continue;
    sw += wi; swx += wi * x[i]; swy += wi * y[i];
    swxx += wi * x[i] * x[i]; swxy += wi * x[i] * y[i];
  }
  double den = sw * swxx - swx * swx;
  double b = (sw * swxy - swx * swy) / den;
  double a = (swy - b * swx) / sw;
  return List::create(_["intercept"] = a, _["slope"] = b);
}

//' Pairwise correlations of matrix columns
//'
//' @param m dense numeric matrix
//' @return correlation matrix
//' @examples
//' m <- matrix(rnorm(30), ncol = 3)
//' lv_col_cor(m)
//' @export
// [[Rcpp::export]]
arma::mat lv_col_cor(const arma::mat& m) {
  return arma::cor(m, 1);
}
