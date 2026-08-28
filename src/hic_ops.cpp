// Insulation-score TAD boundary metrics and boundary-window contact counting.
#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp;

//' Diamond insulation score
//'
//' Mean contact frequency in the lower diamond of radius \code{window}
//' around each bin of a contact matrix; bins near the matrix edges are NA.
//' Local minima of this score mark candidate TAD boundaries.
//'
//' @param m square dense contact matrix
//' @param window diamond radius in bins
//' @return numeric vector of insulation scores per bin
//' @examples
//' m <- matrix(1, 10, 10)            # uniform contacts: flat profile
//' lv_insulation_score(m, window = 2)[3:8]
//' @export
// [[Rcpp::export]]
NumericVector lv_insulation_score(const NumericMatrix& m, int window) {
  int n = m.nrow();
  NumericVector out(n, NA_REAL);
  for (int i = window; i < n - window; i++) {
    double s = 0; int cnt = 0;
    for (int a = i - window; a < i; a++)
      for (int b = i + 1; b <= i + window; b++) { s += m(a, b); cnt++; }
    out[i] = cnt ? s / cnt : NA_REAL;
  }
  return out;
}

//' Total contacts with at least one endpoint inside a boundary window
//'
//' The window spans [bin_start - flank, bin_start + bin_size + flank);
//' bins are 1-based, positions derived as (bin - 1) * bin_size.
//'
//' @param bin1,bin2 integer vectors of pixel endpoint bins
//' @param count pixel counts
//' @param bin_start genomic start (bp) of the boundary bin
//' @param bin_size bin resolution (bp)
//' @param flank window extension on both sides (bp)
//' @return total contact count touching the window
//' @examples
//' lv_window_contacts(c(1L, 4L, 9L), c(2L, 5L, 10L), c(2, 3, 4),
//'                    bin_start = 4000, bin_size = 1000, flank = 1500)
//' @export
// [[Rcpp::export]]
double lv_window_contacts(IntegerVector bin1, IntegerVector bin2,
                          NumericVector count, double bin_start,
                          double bin_size, double flank) {
  double lo = bin_start - flank, hi = bin_start + bin_size + flank;
  double s = 0;
  for (int i = 0; i < bin1.size(); i++) {
    double p1 = (bin1[i] - 1) * bin_size, p2 = (bin2[i] - 1) * bin_size;
    if ((p1 >= lo && p1 < hi) || (p2 >= lo && p2 < hi)) s += count[i];
  }
  return s;
}
