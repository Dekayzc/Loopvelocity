test_that("project_cells applies Eq 16 with EM parameters", {
  d <- make_synthetic_two_sided()
  vel <- list(I = d$I, S = d$S, current = d$I,
              mode = setNames("two-sided", "tad1"),
              par = matrix(c(0.7, 0.25), nrow = 1,
                           dimnames = list("tad1", c("theta", "gamma"))))
  p <- project_cells(vel, delta_t = 1)
  expect_equal(dim(p$projected), dim(d$I))
  expect_equal(as.numeric(p$projected[1, ]),
               as.numeric(d$I[1, ]) +
                 (1 / 0.25 - as.numeric(d$I[1, ])) * (1 - exp(-0.3 * 0.25 * 1)),
               tolerance = 1e-8)
})

test_that("one-sided EM projection uses Eq 6 forward evaluation", {
  t <- runif(10, 0.5, 8)
  sol <- lv_one_sided_solution(t, 1, 0.3, 0.5)
  vel <- list(I = matrix(sol$I, 1, dimnames = list("tad1", paste0("c", 1:10))),
              S = matrix(sol$S, 1, dimnames = list("tad1", paste0("c", 1:10))),
              current = matrix(sol$I, 1, dimnames = list("tad1", paste0("c", 1:10))),
              mode = setNames("one-sided", "tad1"),
              par = matrix(c(0.3, 0.5), nrow = 1,
                           dimnames = list("tad1", c("beta", "gamma"))),
              time = matrix(t, 1, dimnames = list("tad1", paste0("c", 1:10))))
  p <- project_cells(vel, delta_t = 2)
  fwd <- lv_one_sided_solution(t + 2, 1, 0.3, 0.5)$I
  expect_equal(as.numeric(p$projected[1, ]), fwd, tolerance = 1e-8)
})

test_that("steady-state result passes through unchanged", {
  d <- make_synthetic_two_sided(n_cells = 20)
  vel <- loop_velocity_estimates(d$S, d$I, mode = "two-sided", kCells = 1,
                                 verbose = FALSE)
  p <- project_cells(vel)
  expect_equal(p$projected, vel$projected)
})

test_that("projection from EM output is coherent with Eq 12 increment", {
  d <- make_synthetic_two_sided(n_cells = 40, seed = 21)
  vel <- loop_velocity_estimates(d$S, d$I, mode = "two-sided", kCells = 1,
                                 verbose = FALSE)
  vel <- loop_velocity_em(vel, iteration = 10, n.cores = 1, verbose = FALSE)
  p <- project_cells(vel, delta_t = 0.5)
  inc <- as.numeric(p$projected[1, ]) - as.numeric(vel$I[1, ])
  # Eq 16 with EM parameters: increments must match (1/gamma - I)(1 - e^...)
  th <- vel$par["tad1", "theta"]; gm <- vel$par["tad1", "gamma"]
  expect_equal(inc, (1 / gm - as.numeric(vel$I[1, ])) * (1 - exp(-(1 - th) * gm * 0.5)),
               tolerance = 1e-8)
})
