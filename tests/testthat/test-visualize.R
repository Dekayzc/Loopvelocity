test_that("plot_loop_velocity renders both kernels without error", {
  d <- steady_line_two_sided(theta = 0.7, gamma = 0.25, n = 30, n_tads = 3, seed = 1)
  vel <- loop_velocity_estimates(d$S, d$I, mode = "two-sided", kCells = 1,
                                 verbose = FALSE)
  emb <- cbind(runif(30), runif(30)); rownames(emb) <- colnames(d$S)
  cc <- setNames(rep("grey50", 30), colnames(d$S))
  tmp <- tempfile(fileext = ".png"); grDevices::png(tmp)
  expect_error(plot_loop_velocity(emb, vel, kernel = "correlation", n = 10,
                                  cell.colors = cc), NA)
  expect_error(plot_loop_velocity(emb, vel, kernel = "euclidean", n = 10,
                                  cell.colors = cc), NA)
  grDevices::dev.off()
})

test_that("grid flow returns arrow coordinates", {
  d <- steady_line_two_sided(theta = 0.7, gamma = 0.25, n = 30, n_tads = 3, seed = 1)
  vel <- loop_velocity_estimates(d$S, d$I, mode = "two-sided", kCells = 1,
                                 verbose = FALSE)
  emb <- cbind(runif(30), runif(30)); rownames(emb) <- colnames(d$S)
  cc <- setNames(rep("grey50", 30), colnames(d$S))
  tmp <- tempfile(fileext = ".png"); grDevices::png(tmp)
  det <- plot_loop_velocity(emb, vel, kernel = "correlation", n = 10,
                            cell.colors = cc, show.grid.flow = TRUE, grid.n = 5)
  grDevices::dev.off()
  expect_true(is.matrix(det$garrows))
  expect_true(all(c("x0", "y0", "x1", "y1") %in% colnames(det$garrows)))
  expect_equal(dim(det$tp), c(30, 30))
})

test_that("filter_tads_by_activity keeps active TADs", {
  m <- matrix(c(10, 0, 0, 0), nrow = 2, dimnames = list(c("t1", "t2"), c("a", "b")))
  keep <- filter_tads_by_activity(list(S = m), clusters = c(a = "x", b = "y"))
  expect_equal(keep, "t1")
})

test_that("plot_velocity_pcs draws without error", {
  d <- steady_line_two_sided(theta = 0.7, gamma = 0.25, n = 25, n_tads = 3, seed = 4)
  vel <- loop_velocity_estimates(d$S, d$I, mode = "two-sided", kCells = 1,
                                 verbose = FALSE)
  tmp <- tempfile(fileext = ".png"); grDevices::png(tmp)
  cc <- setNames(rep("grey50", 25), colnames(d$S))
  expect_error(plot_velocity_pcs(vel, nPcs = 3, cell.colors = cc), NA)
  grDevices::dev.off()
})

test_that("plot_si_fit reports fitted coefficients", {
  d <- steady_line_two_sided(theta = 0.7, gamma = 0.25, n = 25, n_tads = 30,
                             seed = 3)
  vel <- loop_velocity_estimates(d$S, d$I, mode = "two-sided", kCells = 1,
                                 mult = 1e9, verbose = FALSE)
  tmp <- tempfile(fileext = ".png"); grDevices::png(tmp)
  cf <- plot_si_fit(vel, "tad1")
  grDevices::dev.off()
  expect_equal(cf[["slope"]], (1 - 0.7) / 0.7, tolerance = 0.1)
})
