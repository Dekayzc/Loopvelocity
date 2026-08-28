# Cell-state extrapolation (supplementary notes, Eq 16 two-sided; Eq 5-6
# evaluated forward for one-sided EM fits).

#' Project cell chromatin states forward
#'
#' Two-sided (steady state or EM): dI = (1/gamma - I) * (1 - exp(-(1-theta)
#' * gamma * dt)) with alpha == 1 (Eq 16). One-sided EM: evaluates the
#' analytic solution I(t_c + dt) per cell (Eq 6) at each cell's hidden
#' time. Steady-state results already carry a projection (Eq 16 / Eq 5-6
#' applied during fitting) and pass through unchanged.
#'
#' @param vel velocity result from [loop_velocity_estimates] and optionally
#'   [loop_velocity_em]
#' @param delta_t extrapolation horizon (arbitrary time units)
#' @return the input list with (updated) \code{$projected}
#' @examples
#' set.seed(1); n <- 30
#' t <- sort(runif(n, 0.5, 8))
#' sol <- lv_two_sided_solution(t, 1, 0.7, 0.25)
#' vel <- list(I = matrix(sol$I, 1, dimnames = list("tad1", paste0("c", 1:n))),
#'             S = matrix(sol$S, 1, dimnames = list("tad1", paste0("c", 1:n))),
#'             mode = setNames("two-sided", "tad1"),
#'             par = matrix(c(0.7, 0.25), nrow = 1,
#'                          dimnames = list("tad1", c("theta", "gamma"))))
#' p <- project_cells(vel, delta_t = 1)
#' range(p$projected - vel$I)   # Eq 16 increments
#' @export
project_cells <- function(vel, delta_t = 1) {
  if (!is.null(vel$par)) {   # EM fit present (possibly on a TAD subset)
    has_em <- rownames(vel$I) %in% rownames(vel$par)
    inc <- t(vapply(rownames(vel$I), function(k) {
      if (!k %in% rownames(vel$par)) {
        # not EM-refined: keep the steady-state projection for this TAD
        return(vel$projected[k, ] - vel$I[k, ])
      }
      if (vel$mode[k] == "two-sided") {
        # Eq 16 increment, added to the current state
        as.numeric(lv_project_two_sided(as.numeric(vel$I[k, ]),
                                        vel$par[k, "theta"],
                                        vel$par[k, "gamma"], delta_t, alpha = 1))
      } else {
        # Eq 5-6 forward evaluation: the projected level itself
        as.numeric(lv_one_sided_solution(vel$time[k, ] + delta_t, 1,
                                         vel$par[k, "beta"],
                                         vel$par[k, "gamma"])$I) -
          as.numeric(vel$I[k, ])
      }
    }, numeric(ncol(vel$I))))
    rownames(inc) <- rownames(vel$I); colnames(inc) <- colnames(vel$I)
    proj <- vel$I + inc; proj[proj < 0] <- 0
    vel$projected <- proj
    vel$delta_t <- delta_t
  }
  vel
}
