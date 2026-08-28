// Loop velocity kinetic models: analytic solutions of the loop extrusion
// ODEs and EM objectives (supplementary notes Eq 1-12, 16-18).
#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp;

// Eq 11-12: two-sided loop extrusion solutions with I0 = S0 = 0.
// dS/dt = theta*alpha - theta*gamma*S;  dI/dt = (1-theta)*alpha - (1-theta)*gamma*I
// [[Rcpp::export]]
List lv_two_sided_solution(NumericVector t, double alpha, double theta, double gamma) {
  NumericVector S(t.size()), I(t.size());
  for (int i = 0; i < t.size(); i++) {
    S[i] = alpha / gamma * (1.0 - std::exp(-theta * gamma * t[i]));
    I[i] = alpha / gamma * (1.0 - std::exp(-(1.0 - theta) * gamma * t[i]));
  }
  List out = List::create(_["S"] = S, _["I"] = I);
  return out;
}

// Eq 5-6: one-sided loop extrusion solutions with I0 = S0 = 0.
// dS/dt = alpha - beta*S;  dI/dt = beta*S - gamma*I
// The alpha/(gamma-beta)*(e^{-gamma t} - e^{-beta t}) term of Eq 6 tends to
// -alpha * t * e^{-beta t} as gamma -> beta (l'Hopital limit); guarded to
// avoid catastrophic cancellation near the degenerate case.
// [[Rcpp::export]]
List lv_one_sided_solution(NumericVector t, double alpha, double beta, double gamma) {
  NumericVector S(t.size()), I(t.size());
  for (int i = 0; i < t.size(); i++) {
    S[i] = alpha / beta * (1.0 - std::exp(-beta * t[i]));
    double base = alpha / gamma * (1.0 - std::exp(-gamma * t[i]));
    if (std::fabs(gamma - beta) < 1e-6 * (gamma + beta)) {
      I[i] = base - alpha * t[i] * std::exp(-beta * t[i]);
    } else {
      I[i] = base + alpha / (gamma - beta) *
        (std::exp(-gamma * t[i]) - std::exp(-beta * t[i]));
    }
  }
  List out = List::create(_["S"] = S, _["I"] = I);
  return out;
}

// Eq 18 negative log-likelihood shape (up to constants), two-sided
// parameterization. phi = c(alpha, theta, gamma).
// [[Rcpp::export]]
double lv_em_objective_two_sided(NumericVector t, NumericVector Sobs,
                                 NumericVector Iobs, NumericVector phi) {
  List sol = lv_two_sided_solution(t, phi[0], phi[1], phi[2]);
  NumericVector Sm = sol["S"], Im = sol["I"];
  double ss = 0.0;
  for (int i = 0; i < t.size(); i++)
    ss += (Iobs[i] - Im[i]) * (Iobs[i] - Im[i]) + (Sobs[i] - Sm[i]) * (Sobs[i] - Sm[i]);
  return ss;
}

// Eq 18 objective, one-sided parameterization. phi = c(alpha, beta, gamma).
// [[Rcpp::export]]
double lv_em_objective_one_sided(NumericVector t, NumericVector Sobs,
                                 NumericVector Iobs, NumericVector phi) {
  List sol = lv_one_sided_solution(t, phi[0], phi[1], phi[2]);
  NumericVector Sm = sol["S"], Im = sol["I"];
  double ss = 0.0;
  for (int i = 0; i < t.size(); i++)
    ss += (Iobs[i] - Im[i]) * (Iobs[i] - Im[i]) + (Sobs[i] - Sm[i]) * (Sobs[i] - Sm[i]);
  return ss;
}

// Eq 16: two-sided projection horizon delta over observed I values.
// dI = (alpha/gamma - I) * (1 - exp(-(1-theta)*gamma*delta_t)).
// [[Rcpp::export]]
NumericVector lv_project_two_sided(NumericVector Iobs, double theta, double gamma,
                                   double delta_t, double alpha) {
  NumericVector out(Iobs.size());
  for (int i = 0; i < Iobs.size(); i++)
    out[i] = (alpha / gamma - Iobs[i]) * (1.0 - std::exp(-(1.0 - theta) * gamma * delta_t));
  return out;
}
