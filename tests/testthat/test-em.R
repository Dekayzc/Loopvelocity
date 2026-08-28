test_that("EM recovers two-sided phi from synthetic trajectories", {
  d <- make_synthetic_two_sided(n_cells = 80, alpha = 1, theta = 0.7,
                                gamma = 0.25, noise = 0.02, seed = 7)
  vel <- list(S = d$S, I = d$I,
              ko = data.frame(theta = d$truth[["theta"]], row.names = "tad1",
                              gamma = d$truth[["gamma"]]),
              mode = setNames("two-sided", "tad1"))
  em <- loop_velocity_em(vel, iteration = 40, tol = 0.01, n.cores = 1)
  expect_equal(em$par["tad1", "theta"], 0.7, tolerance = 0.08)
  expect_equal(em$par["tad1", "gamma"], 0.25, tolerance = 0.08)
})

test_that("EM converges before max iterations on clean data", {
  d <- make_synthetic_two_sided(n_cells = 80, noise = 0.005, seed = 8)
  vel <- list(S = d$S, I = d$I,
              ko = data.frame(theta = d$truth[["theta"]], row.names = "tad1",
                              gamma = d$truth[["gamma"]]),
              mode = setNames("two-sided", "tad1"))
  em <- loop_velocity_em(vel, iteration = 60, tol = 0.01, n.cores = 1)
  expect_true(em$converged[["tad1"]])
  expect_lt(em$iterations[["tad1"]], 60)
})

test_that("EM recovers one-sided phi from synthetic trajectories", {
  # alpha = 1: the EM fixes the rate scale with unit TAD generation rate
  d <- make_synthetic_one_sided(n_cells = 80, alpha = 1, beta = 0.3,
                                gamma = 0.5, noise = 0.01, seed = 9)
  vel <- list(S = d$S, I = d$I,
              ko = data.frame(theta = NA_real_, gamma = NA_real_, row.names = "tad1",
                              ratio = d$truth[["gamma"]] / d$truth[["beta"]]),
              mode = setNames("one-sided", "tad1"))
  em <- loop_velocity_em(vel, iteration = 40, tol = 0.01, n.cores = 1)
  expect_equal(em$par["tad1", "beta"], 0.3, tolerance = 0.1)
  expect_equal(em$par["tad1", "gamma"], 0.5, tolerance = 0.1)
})

test_that("hidden times increase with observed I", {
  d <- make_synthetic_two_sided(n_cells = 80, noise = 0.01, seed = 10)
  vel <- list(S = d$S, I = d$I,
              ko = data.frame(theta = d$truth[["theta"]], row.names = "tad1",
                              gamma = d$truth[["gamma"]]),
              mode = setNames("two-sided", "tad1"))
  em <- loop_velocity_em(vel, iteration = 20, tol = 0.01, n.cores = 1)
  expect_gt(cor(em$time["tad1", ], d$t), 0.9)
})
