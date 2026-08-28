# EM dynamical model (supplementary notes, "Dynamical model", Eq 17-19).
# Hidden variable: per-cell time t. E-step: t_c = argmin over t of
#   (I_c - I(t, phi))^2 + (S_c - S(t, phi))^2   (Eq 19)
# via box-constrained L-BFGS. M-step: the same objective minimized over
# phi, again with L-BFGS. Iterated until the maximum relative parameter
# change falls below tol or the iteration cap is hit (notes: 100
# iterations or < 1% change).

# invert I(t) (Eq 12) for an initial time guess, two-sided
.init_time_two_sided <- function(Iobs, alpha, theta, gamma) {
  sat <- alpha / gamma
  frac <- pmin(pmax(Iobs / sat, 0), 1 - 1e-6)
  -log(1 - frac) / ((1 - theta) * gamma)
}

# invert I(t) (Eq 6) numerically, one-sided
.init_time_one_sided <- function(Iobs, alpha, beta, gamma) {
  vapply(Iobs, function(iv) {
    stats::optimize(function(tt)
      (lv_one_sided_solution(c(tt), alpha, beta, gamma)$I[1] - iv)^2,
      interval = c(0, 50 / min(beta, gamma)))$minimum
  }, numeric(1))
}

# The absolute rate scale is unidentifiable when times are hidden:
# (alpha, beta, gamma, t) -> (c*alpha, c*beta, c*gamma, t/c) yields
# identical trajectories. Following the notes' convention ("Let alpha = 1,
# the generation rate of TADs in DNA is unit rate 1"), alpha is fixed at 1
# and the fitted phi is c(theta, gamma) or c(beta, gamma).
.em_fit_one_tad <- function(s, i, init, mode, iteration, tol) {
  two_sided <- mode == "two-sided"
  par_box <- if (two_sided)
    list(lower = c(0.51, 1e-6), upper = c(0.99, 1e2))
  else list(lower = c(1e-6, 1e-6), upper = c(1e2, 1e2))
  obj <- if (two_sided)
    function(par, tt) lv_em_objective_two_sided(tt, s, i, c(1, par[1], par[2]))
  else function(par, tt) lv_em_objective_one_sided(tt, s, i, c(1, par[1], par[2]))

  par <- init
  tt <- if (two_sided)
    .init_time_two_sided(i, 1, par[1], par[2])
  else .init_time_one_sided(i, 1, par[1], par[2])

  rel <- Inf; it <- 0L
  while (it < iteration && rel > tol) {
    it <- it + 1L
    # M-step: parameters from current hidden times
    om <- stats::optim(par, function(p) obj(p, tt), method = "L-BFGS-B",
                       lower = par_box$lower, upper = par_box$upper)
    par_new <- om$par
    # E-step: per-cell hidden times from current parameters (Eq 19) --
    # each cell minimizes its own residual, not the pooled one
    tmax <- if (two_sided) 20 / ((1 - par_new[1]) * par_new[2])
            else 20 / min(par_new[1], par_new[2])
    tt <- vapply(seq_along(s), function(k) {
      stats::optim(max(tt[k], 1e-3), function(tk) {
        sol <- if (two_sided) lv_two_sided_solution(c(tk), 1, par_new[1], par_new[2])
               else lv_one_sided_solution(c(tk), 1, par_new[1], par_new[2])
        (s[k] - sol$S[1])^2 + (i[k] - sol$I[1])^2
      }, method = "L-BFGS-B", lower = 0, upper = max(tmax, 1))$par
    }, numeric(1))
    rel <- max(abs(par_new - par) / pmax(abs(par), 1e-8))
    par <- par_new
  }
  names(par) <- if (two_sided) c("theta", "gamma") else c("beta", "gamma")
  list(par = par, t = tt, iterations = it, converged = rel <= tol)
}

#' EM dynamical-model loop velocity
#'
#' Refines steady-state parameters with an EM algorithm over the full
#' loop-extrusion ODE solutions (Eq 11-12 / Eq 5-6), treating per-cell
#' time as the hidden variable (Eq 17-19). Both the E-step (hidden times)
#' and the M-step (parameters) use box-constrained L-BFGS, as specified
#' in the supplementary notes.
#'
#' @param vel output of [loop_velocity_estimates]: its normalized \code{$S},
#'   \code{$I}, \code{$mode} and \code{$ko} fields seed the EM
#' @param iteration iteration cap (notes default 100)
#' @param tol convergence: maximum relative parameter change between rounds
#' @param n.cores worker processes (per-TAD parallelism)
#' @param verbose progress messages
#' @return the input list augmented with \code{$par} (TAD x phi), \code{$time}
#'   (TAD x cell hidden times), \code{$iterations}, \code{$converged}
#' @examples
#' # simulate cells along the two-sided trajectory and recover theta/gamma
#' set.seed(1); n <- 40
#' t <- sort(runif(n, 0.5, 10))
#' th <- 0.7; gm <- 0.3
#' sol <- lv_two_sided_solution(t, alpha = 1, theta = th, gamma = gm)
#' vel <- list(S = matrix(sol$S, 1, dimnames = list("tad1", paste0("c", 1:n))),
#'             I = matrix(sol$I, 1, dimnames = list("tad1", paste0("c", 1:n))),
#'             mode = setNames("two-sided", "tad1"),
#'             ko = data.frame(theta = th, gamma = gm, ratio = NA,
#'                             geff = (1 - th) / th, o = 0, row.names = "tad1"))
#' em <- loop_velocity_em(vel, iteration = 20, n.cores = 1)
#' em$par             # close to theta = 0.7, gamma = 0.3
#' @export
loop_velocity_em <- function(vel, iteration = 100, tol = 0.01,
                             n.cores = default_n_cores(), verbose = TRUE) {
  vt <- rownames(vel$S)
  if (verbose) cat("EM dynamical model on", length(vt), "TADs ... ")
  res <- parallel::mclapply(vt, function(k) {
    if (vel$mode[k] == "two-sided") {
      th0 <- suppressWarnings(vel$ko[k, "theta"]); gm0 <- suppressWarnings(vel$ko[k, "gamma"])
      if (is.na(th0)) th0 <- 0.6
      if (is.na(gm0) || gm0 <= 0) gm0 <- 0.1
      init <- c(theta = min(max(th0, 0.55), 0.95), gamma = min(max(gm0, 1e-3), 1e2))
    } else {
      rt0 <- suppressWarnings(vel$ko[k, "ratio"])
      if (is.na(rt0) || rt0 <= 0) rt0 <- 1
      # data-driven init: S saturates at alpha/beta, so beta0 = 1/max(S);
      # the steady-state ratio then fixes gamma0 = (gamma/beta) * beta0
      smax <- max(vel$S[k, ])
      b0 <- if (is.finite(smax) && smax > 0) 1 / smax else 1
      init <- c(beta = min(max(b0, 1e-3), 1e2),
                gamma = min(max(rt0 * b0, 1e-3), 1e2))
    }
    .em_fit_one_tad(as.numeric(vel$S[k, ]), as.numeric(vel$I[k, ]),
                    init, vel$mode[k], iteration, tol)
  }, mc.cores = n.cores, mc.preschedule = TRUE)
  names(res) <- vt
  par <- do.call(rbind, lapply(res, `[[`, "par"))
  time <- do.call(rbind, lapply(res, `[[`, "t"))
  rownames(par) <- rownames(time) <- vt
  colnames(time) <- colnames(vel$S)
  vel$par <- par; vel$time <- time
  vel$iterations <- vapply(res, function(x) x$iterations, numeric(1))
  vel$converged <- vapply(res, function(x) x$converged, logical(1))
  if (verbose) cat("done\n")
  vel
}
