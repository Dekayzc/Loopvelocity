# Shared synthetic-data generator: draws cell times, evaluates the
# analytic ODE solutions (ground truth), adds multiplicative noise.

make_synthetic_two_sided <- function(n_cells = 60, alpha = 1, theta = 0.7,
                                     gamma = 0.25, tmax = 12, noise = 0.01,
                                     seed = 1, tad_name = "tad1") {
  set.seed(seed)
  t <- sort(stats::runif(n_cells, 0.2, tmax))
  sol <- lv_two_sided_solution(t, alpha, theta, gamma)
  cn <- paste0("c", seq_len(n_cells))
  list(t = t,
       S = matrix(pmax(0, sol$S * (1 + stats::rnorm(n_cells, 0, noise))), nrow = 1,
                  dimnames = list(tad_name, cn)),
       I = matrix(pmax(0, sol$I * (1 + stats::rnorm(n_cells, 0, noise))), nrow = 1,
                  dimnames = list(tad_name, cn)),
       truth = c(alpha = alpha, theta = theta, gamma = gamma))
}

make_synthetic_one_sided <- function(n_cells = 60, alpha = 1, beta = 0.2,
                                     gamma = 0.4, tmax = 12, noise = 0.01,
                                     seed = 1, tad_name = "tad1") {
  set.seed(seed)
  t <- sort(stats::runif(n_cells, 0.2, tmax))
  sol <- lv_one_sided_solution(t, alpha, beta, gamma)
  cn <- paste0("c", seq_len(n_cells))
  list(t = t,
       S = matrix(pmax(0, sol$S * (1 + stats::rnorm(n_cells, 0, noise))), nrow = 1,
                  dimnames = list(tad_name, cn)),
       I = matrix(pmax(0, sol$I * (1 + stats::rnorm(n_cells, 0, noise))), nrow = 1,
                  dimnames = list(tad_name, cn)),
       truth = c(alpha = alpha, beta = beta, gamma = gamma))
}

steady_line_two_sided <- function(theta, gamma, n = 60, n_tads = 4, noise = 0.01,
                                  seed = 1) {
  set.seed(seed)
  a <- (1 - theta) / theta
  b <- (1 / gamma) * (2 - 1 / theta)
  cn <- paste0("c", seq_len(n)); tn <- paste0("tad", seq_len(n_tads))
  per_tad <- lapply(tn, function(g) {
    Ig <- runif(n, 0.5, 3)
    list(I = Ig, S = pmax(0, a * Ig + b + rnorm(n, 0, noise * b)))
  })
  S <- t(sapply(per_tad, `[[`, "S")); I <- t(sapply(per_tad, `[[`, "I"))
  dimnames(S) <- dimnames(I) <- list(tn, cn)
  list(S = S, I = I, truth = c(theta = theta, gamma = gamma))
}

steady_line_one_sided <- function(ratio, n = 60, n_tads = 4, noise = 0.01, seed = 1) {
  set.seed(seed)
  cn <- paste0("c", seq_len(n)); tn <- paste0("tad", seq_len(n_tads))
  per_tad <- lapply(tn, function(g) {
    Ig <- runif(n, 0.5, 3)
    list(I = Ig, S = pmax(0, ratio * Ig * (1 + rnorm(n, 0, noise))))
  })
  S <- t(sapply(per_tad, `[[`, "S")); I <- t(sapply(per_tad, `[[`, "I"))
  dimnames(S) <- dimnames(I) <- list(tn, cn)
  list(S = S, I = I, truth = c(ratio = ratio))
}
