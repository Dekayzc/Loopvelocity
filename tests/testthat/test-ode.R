sim_two_sided <- function(t, alpha, theta, gamma) {
  list(S = alpha / gamma * (1 - exp(-theta * gamma * t)),
       I = alpha / gamma * (1 - exp(-(1 - theta) * gamma * t)))
}
sim_one_sided <- function(t, alpha, beta, gamma) {
  S <- alpha / beta * (1 - exp(-beta * t))
  I <- alpha / gamma * (1 - exp(-gamma * t)) +
       alpha / (gamma - beta) * (exp(-gamma * t) - exp(-beta * t))
  list(S = S, I = I)
}

test_that("two-sided ODE solutions match Eq 11-12", {
  tt <- seq(0, 50, length.out = 7)
  got <- lv_two_sided_solution(tt, alpha = 2, theta = 0.7, gamma = 0.3)
  ref <- sim_two_sided(tt, 2, 0.7, 0.3)
  expect_equal(got$S, ref$S, tolerance = 1e-8)
  expect_equal(got$I, ref$I, tolerance = 1e-8)
  expect_equal(unlist(lv_two_sided_solution(0, 2, 0.7, 0.3)), c(S = 0, I = 0))
})

test_that("one-sided ODE solutions match Eq 5-6 incl. gamma~beta limit", {
  tt <- seq(0, 40, length.out = 9)
  got <- lv_one_sided_solution(tt, alpha = 1.5, beta = 0.2, gamma = 0.4)
  ref <- sim_one_sided(tt, 1.5, 0.2, 0.4)
  expect_equal(got$S, ref$S, tolerance = 1e-8)
  expect_equal(got$I, ref$I, tolerance = 1e-8)
  # degenerate gamma == beta: series limit of Eq 6
  gb <- lv_one_sided_solution(tt, alpha = 1.5, beta = 0.3, gamma = 0.3 + 1e-4)
  gd <- lv_one_sided_solution(tt, alpha = 1.5, beta = 0.3, gamma = 0.3)
  expect_equal(gb$I, gd$I, tolerance = 1e-3)
})

test_that("EM objective Eq 18 sums squared residuals over cells", {
  Sobs <- c(1, 2, 3); Iobs <- c(0.5, 1.1, 1.8); tt <- c(1, 2, 3)
  phi <- c(alpha = 1, theta = 0.6, gamma = 0.2)
  expect_equal(
    lv_em_objective_two_sided(tt, Sobs, Iobs, phi),
    sum((Iobs - sim_two_sided(tt, 1, 0.6, 0.2)$I)^2 +
        (Sobs - sim_two_sided(tt, 1, 0.6, 0.2)$S)^2),
    tolerance = 1e-10)
})

test_that("one-sided EM objective matches R reference", {
  Sobs <- c(1, 2, 3); Iobs <- c(0.5, 1.1, 1.8); tt <- c(1, 2, 3)
  phi <- c(1.2, 0.25, 0.5)
  ref <- sim_one_sided(tt, 1.2, 0.25, 0.5)
  expect_equal(
    lv_em_objective_one_sided(tt, Sobs, Iobs, phi),
    sum((Iobs - ref$I)^2 + (Sobs - ref$S)^2),
    tolerance = 1e-10)
})

test_that("two-sided projection delta matches Eq 16", {
  Iobs <- c(0.2, 0.8, 2.0)
  got <- lv_project_two_sided(Iobs, theta = 0.6, gamma = 0.3, delta_t = 2, alpha = 1)
  ref <- (1 / 0.3 - Iobs) * (1 - exp(-(1 - 0.6) * 0.3 * 2))
  expect_equal(got, ref, tolerance = 1e-12)
})
