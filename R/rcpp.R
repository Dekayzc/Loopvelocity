# Package-level import directives and compiled-code registration.
# C++ sources: loop_models.cpp (ODE solutions, EM objectives, projection),
# fast_stats.cpp (weighted fits, correlations), hic_ops.cpp (insulation,
# boundary-window contacts).
NULL

#' Two-sided loop extrusion ODE solutions (Eq 11-12)
#'
#' Evaluates the analytic solutions of dS/dt = theta*alpha - theta*gamma*S
#' and dI/dt = (1-theta)*alpha - (1-theta)*gamma*I with I0 = S0 = 0.
#' @param t numeric vector of time points
#' @param alpha contact formation rate
#' @param theta fraction of the extrusion rate attributed to the S boundary
#' @param gamma loop dissociation rate
#' @return list with elements \code{S} and \code{I}
#' @examples
#' sol <- lv_two_sided_solution(seq(0, 10, by = 2), alpha = 1, theta = 0.7,
#'                              gamma = 0.3)
#' plot(sol$I, sol$S, type = "b")
#' @export
#' @name lv_two_sided_solution
NULL

#' One-sided loop extrusion ODE solutions (Eq 5-6)
#'
#' Evaluates the analytic solutions of dS/dt = alpha - beta*S and
#' dI/dt = beta*S - gamma*I with I0 = S0 = 0; the gamma == beta limit is
#' handled analytically.
#' @param t numeric vector of time points
#' @param alpha contact formation rate
#' @param beta active contact reduction rate at the S boundary
#' @param gamma loop dissociation rate at the I boundary
#' @return list with elements \code{S} and \code{I}
#' @examples
#' sol <- lv_one_sided_solution(seq(0, 10, by = 2), alpha = 1, beta = 0.3,
#'                              gamma = 0.5)
#' plot(sol$I, sol$S, type = "b")
#' @export
#' @name lv_one_sided_solution
NULL

#' Two-sided EM objective (Eq 18)
#'
#' Sums squared residuals between observed (S, I) pairs and the model
#' solutions at times \code{t} for parameters \code{phi}.
#' @param t hidden cell times
#' @param Sobs,Iobs observed boundary contacts
#' @param phi numeric vector c(alpha, theta, gamma)
#' @return sum of squared residuals
#' @examples
#' t <- c(1, 2, 3)
#' lv_em_objective_two_sided(t, Sobs = c(0.2, 0.3, 0.5),
#'                           Iobs = c(0.1, 0.2, 0.4),
#'                           phi = c(1, 0.7, 0.3))
#' @export
#' @name lv_em_objective_two_sided
NULL

#' One-sided EM objective (Eq 18)
#'
#' @param t hidden cell times
#' @param Sobs,Iobs observed boundary contacts
#' @param phi numeric vector c(alpha, beta, gamma)
#' @return sum of squared residuals
#' @examples
#' t <- c(1, 2, 3)
#' lv_em_objective_one_sided(t, Sobs = c(0.2, 0.3, 0.5),
#'                           Iobs = c(0.1, 0.2, 0.4),
#'                           phi = c(1, 0.3, 0.5))
#' @export
#' @name lv_em_objective_one_sided
NULL

#' Two-sided projection horizon (Eq 16)
#'
#' dI = (alpha/gamma - I) * (1 - exp(-(1-theta)*gamma*delta_t)).
#' @param Iobs observed I-boundary contacts
#' @param theta,gamma,alpha model parameters
#' @param delta_t extrapolation horizon
#' @return numeric vector of projected increments
#' @examples
#' lv_project_two_sided(Iobs = c(0.2, 0.8, 2), theta = 0.6, gamma = 0.3,
#'                      delta_t = 2, alpha = 1)
#' @export
#' @name lv_project_two_sided
NULL

#' @useDynLib Loopvelocity, .registration = TRUE
#' @importFrom Rcpp evalCpp
#' @import Matrix
#' @importFrom stats lm optim quantile cor sd setNames p.adjust runif rnorm
#' @importFrom graphics arrows legend lines par plot points polygon
#' @importFrom grDevices adjustcolor colorRampPalette rgb
#' @importFrom utils read.delim
#' @importFrom methods as
NULL
