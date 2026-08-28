test_that("weighted LSQ matches lm on quantile-selected points", {
  set.seed(1)
  x <- runif(200); y <- 0.3 * x + 0.7 + rnorm(200, sd = 0.02)
  w <- as.numeric(x >= quantile(x, 0.95) | x <= quantile(x, 0.05))
  fl <- lv_wls_fit(x, y, w)
  ref <- lm(y ~ x, weights = w)
  expect_equal(fl$intercept, unname(coef(ref)[1]), tolerance = 1e-6)
  expect_equal(fl$slope, unname(coef(ref)[2]), tolerance = 1e-6)
})

test_that("weighted LSQ handles zero weights", {
  set.seed(11)
  x <- runif(50); y <- 1.2 * x + 0.4 + rnorm(50, sd = 0.01)
  w <- rep(1:0, length.out = 50)
  fl <- lv_wls_fit(x, y, w)
  ref <- lm(y ~ x, weights = w)
  expect_equal(fl$slope, unname(coef(ref)[2]), tolerance = 1e-6)
  expect_equal(fl$intercept, unname(coef(ref)[1]), tolerance = 1e-6)
})

test_that("column correlations match cor()", {
  set.seed(2)
  m <- matrix(rnorm(50 * 4), ncol = 4)
  expect_equal(lv_col_cor(m), cor(m), tolerance = 1e-8)
})
