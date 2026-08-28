// Balanced k-nearest-neighbour graph construction over single cells.
#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
#include <sstream>
using namespace Rcpp;

// Given a cell x cell distance matrix, form a balanced kNN assignment: cells
// that appear too often as neighbours are capped (maxl), and cells are
// processed greedily in order of how many others selected them. Used for
// contact-profile pooling before steady-state fitting.
// [[Rcpp::export]]
arma::sp_mat lv_balanced_knn_graph(const arma::mat& d, int k, int maxl,
                                   bool return_distance_values = false,
                                   int nthreads = 1) {
  arma::uvec l(d.n_cols, arma::fill::zeros);   // how often each cell was picked
  arma::umat dsi(d.n_rows, d.n_cols);          // per-column sort order of d
#pragma omp parallel for shared(l) num_threads(nthreads)
  for (int i = 0; i < d.n_cols; i++) {
    arma::uvec si = arma::sort_index(d.col(i));
    dsi.col(i) = si;
    l.elem(si.subvec(0, k - 1)) += 1;
  }

  arma::uvec lsi = arma::sort_index(l, "descend");  // greedy processing order
  l.zeros();

  arma::uvec rowind(k * d.n_cols);
  arma::vec vals(k * d.n_cols, arma::fill::ones);

  for (int i = 0; i < d.n_cols; i++) {
    int el = lsi[i];
    arma::uvec si = dsi.col(el);
    int p = 0, j = 0;
    for (j = 0; j < (int)dsi.n_rows && p < k; j++) {
      int m = si[j];
      if (el == m) continue;          // skip self
      if (l[m] >= maxl) continue;     // neighbour saturated
      rowind[el * k + p] = m; l[m]++; p++;
      if (return_distance_values) vals[el * k + p] = d(m, el);
    }
    if (j == (int)dsi.n_rows && p < k) {
      rowind[el * k + p] = el; p++;
      while (p < k) {
        std::stringstream es;
        es << "not enough distinct neighbours for cell " << el
           << "; filling with self-links";
        Rf_warning(es.str().c_str());
        rowind[el * k + p] = el; p++;
      }
    }
  }

  arma::uvec colptr(d.n_cols + 1);
  int cc = 0;
  for (arma::uword i = 0; i < colptr.n_elem; i++) { colptr[i] = cc; cc += k; }
  arma::sp_mat knn(rowind, colptr, vals, d.n_rows, d.n_cols);
  return knn;
}
