// Embedding-transition kernels for velocity visualization: delta
// correlations between cells and expected arrow displacements.
#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp;

// Correlation between per-cell contact deltas d[,i] and the row-wise
// differences of e relative to cell i, across all cells (linear scale).
// [[Rcpp::export]]
arma::mat lv_delta_cor(const arma::mat& e, const arma::mat& d, int nthreads = 1) {
  arma::mat rm(e.n_cols, e.n_cols);
#pragma omp parallel for shared(rm) num_threads(nthreads)
  for (int i = 0; i < (int)e.n_cols; i++) {
    arma::mat t(e); t.each_col() -= e.col(i);
    rm.col(i) = arma::cor(t, d.col(i));
  }
  return rm;
}

// Same with signed sqrt transform of the differences.
// [[Rcpp::export]]
arma::mat lv_delta_cor_sqrt(const arma::mat& e, const arma::mat& d, int nthreads = 1) {
  arma::mat rm(e.n_cols, e.n_cols);
#pragma omp parallel for shared(rm) num_threads(nthreads)
  for (int i = 0; i < (int)e.n_cols; i++) {
    arma::mat t(e); t.each_col() -= e.col(i);
    t = sqrt(abs(t)) % sign(t);
    rm.col(i) = arma::cor(t, d.col(i));
  }
  return rm;
}

// Same with signed log10(|.| + pseudocount) transform of the differences.
// [[Rcpp::export]]
arma::mat lv_delta_cor_log10(const arma::mat& e, const arma::mat& d,
                             double pseudocount = 1.0, int nthreads = 1) {
  arma::mat rm(e.n_cols, e.n_cols);
#pragma omp parallel for shared(rm) num_threads(nthreads)
  for (int i = 0; i < (int)e.n_cols; i++) {
    arma::mat t(e); t.each_col() -= e.col(i);
    t = log10(abs(t) + pseudocount) % sign(t);
    rm.col(i) = arma::cor(t, d.col(i));
  }
  return rm;
}

// Euclidean distances between the columns of e and the columns of p
// (cells x cells when e and p are TADs x cells).
// [[Rcpp::export]]
arma::mat lv_col_euclid(const arma::mat& e, const arma::mat& p, int nthreads = 1) {
  arma::mat rm(e.n_cols, p.n_cols);
#pragma omp parallel for shared(rm) num_threads(nthreads)
  for (int i = 0; i < (int)p.n_cols; i++) {
    arma::mat t(e); t.each_col() -= p.col(i);
    t %= t;
    arma::rowvec v = sqrt(sum(t, 0));
    rm.col(i) = trans(v);
  }
  return rm;
}

// Expected arrow displacement on an embedding given cell-cell transition
// probabilities: normalized displacement directions weighted by tp, minus
// the same average over the binarized neighborhood (kNN center).
// [[Rcpp::export]]
arma::mat lv_emb_arrows(const arma::mat& emb, const arma::sp_mat& tp,
                        double arrow_scale = 1.0, int nthreads = 1) {
  arma::mat dm(emb.n_cols, emb.n_rows);
  arma::sp_mat tpb(tp);
  arma::vec tprs(tp.n_cols, arma::fill::zeros);
  for (arma::sp_mat::iterator ci = tpb.begin(); ci != tpb.end(); ++ci)
    tprs[ci.col()]++;
  for (arma::sp_mat::iterator ci = tpb.begin(); ci != tpb.end(); ++ci)
    (*ci) = 1.0 / tprs[ci.col()];

  arma::colvec zv(emb.n_cols, arma::fill::zeros);
  arma::mat temb = trans(emb);
#pragma omp parallel for shared(dm) num_threads(nthreads)
  for (int i = 0; i < (int)emb.n_rows; i++) {
    arma::mat di(temb);
    di.each_col() -= di.col(i);
    di = arma::normalise(di, 2, 0) * arrow_scale;
    di.col(i) = zv;
    arma::vec ds = di * tp.col(i) - di * tpb.col(i);
    dm.col(i) = ds;
  }
  return dm;
}
