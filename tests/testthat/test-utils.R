test_that("balanced_knn returns symmetric-normalized kNN weights", {
  set.seed(3)
  m <- matrix(rnorm(60), nrow = 20)  # cells x features
  rownames(m) <- paste0("c", 1:20)
  knn <- balanced_knn(m, k = 5)
  expect_equal(dim(knn), c(20, 20))
  expect_equal(rownames(knn), rownames(m))
  expect_true(all(abs(knn) <= 1 + 1e-8))
  expect_equal(unname(Matrix::rowSums(knn)), rep(5, 20))
})

test_that("balanced_knn caps neighbour popularity", {
  set.seed(31)
  # one cell extremely similar to all others
  m <- matrix(rnorm(60), nrow = 20)
  m[2, ] <- m[1, ]
  m[3, ] <- m[1, ] * 1.0001
  rownames(m) <- paste0("c", 1:20)
  knn <- balanced_knn(m, k = 3, maxl = 3)
  expect_true(all(Matrix::colSums(knn) <= 3))
})

test_that("contact palette helpers return valid colors", {
  p <- lv_expression_palette()
  expect_length(p, 256)
  cols <- lv_val2col(c(-2, -1, 0, 1, 2))
  expect_length(cols, 5)
})
