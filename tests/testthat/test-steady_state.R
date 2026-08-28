# Steady-state estimator tests. The steady-state model treats the cell
# population as an ensemble on the Eq 13/14 constraint lines, so the
# ground-truth fixtures generate cells ON those lines with noise. (Full
# ODE-trajectory recovery is covered by the EM tests.)
#
# Fixtures use several TADs: with a single TAD, per-cell size
# normalization (colSums) is degenerate -- the column sum is the value
# itself and every normalized cell becomes a constant.


test_that("steady-state fit recovers theta and gamma (Eq 14, alpha=1)", {
  d <- steady_line_two_sided(theta = 0.7, gamma = 0.25, n_tads = 30, seed = 11)
  # mult = average per-cell total contacts keeps cs ~ 1, so gamma stays in
  # the fixture's raw units (alpha == 1 convention ties gamma's scale to
  # the normalization)
  # mult = 1e9 clamps both size factors to the same constant: the
  # normalization degenerates to uniform scaling and the regression
  # parameters are recovered exactly (on real data the per-matrix factors
  # follow the legacy velocity-framework preprocessing)
  fit <- loop_velocity_estimates(d$S, d$I, mode = "two-sided", kCells = 1,
                                 mult = 1e9)
  expect_equal(mean(fit$theta), 0.7, tolerance = 0.05)
  # gamma carries the uniform 1e-3 size-factor scale (alpha = 1 convention)
  expect_equal(mean(fit$gamma, na.rm = TRUE), 0.25 * 1e-3, tolerance = 1.5e-4)
})

test_that("steady-state one-sided recovers gamma/beta ratio (Eq 13)", {
  d <- steady_line_one_sided(ratio = 2, seed = 5)
  fit <- loop_velocity_estimates(d$S, d$I, mode = "one-sided", kCells = 1,
                                 mult = 1e9, fit.quantile = 0.05)
  expect_equal(mean(fit$ratio), 2, tolerance = 0.15)
})

test_that("steady state holds asymptotically for two-sided trajectories", {
  th <- 0.7; gm <- 0.25
  sol <- lv_two_sided_solution(rep(60, 50), 1, th, gm)
  a <- (1 - th) / th; b <- (2 - 1 / th) / gm
  expect_equal(sol$S[1], a * sol$I[1] + b, tolerance = 0.05)
})

test_that("result structure is complete", {
  d <- steady_line_two_sided(theta = 0.7, gamma = 0.25, n = 30, seed = 2)
  fit <- loop_velocity_estimates(d$S, d$I, mode = "two-sided", kCells = 1)
  expect_true(all(c("theta", "gamma", "ratio", "mode", "ko",
                    "current", "projected", "deltaE") %in% names(fit)))
  expect_equal(dim(fit$projected), dim(d$S))
  expect_false(any(fit$projected < 0))
})

test_that("cell pooling via kCells=5 keeps parameter recovery", {
  d <- steady_line_two_sided(theta = 0.7, gamma = 0.25, n = 60, n_tads = 6, seed = 4)
  fit <- loop_velocity_estimates(d$S, d$I, mode = "two-sided", kCells = 5,
                                 mult = 1e9, fit.quantile = 0.05)
  expect_equal(mean(fit$theta), 0.7, tolerance = 0.06)
})

test_that("auto mode uses per-TAD modes from assign_si_boundaries output", {
  d <- steady_line_two_sided(theta = 0.7, gamma = 0.25, seed = 1)
  mtbl <- list(S = d$S, mode = setNames(rep("one-sided", nrow(d$S)), rownames(d$S)))
  class(mtbl) <- "lv_si"
  fit <- loop_velocity_estimates(d$S, d$I, mode = "auto", kCells = 1, mode_tbl = mtbl)
  expect_true(all(as.character(fit$mode) == "one-sided"))
  expect_true(all(!is.na(fit$ratio)))
})
