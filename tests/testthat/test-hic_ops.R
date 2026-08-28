test_that("insulation score dips at a domain boundary", {
  n <- 40; m <- matrix(0, n, n)
  # two domains: strong contacts within bins 1-20 and 21-40, none across
  idx <- expand.grid(i = 1:n, j = 1:n)
  same <- (idx$i <= 20) == (idx$j <= 20)
  m[as.matrix(idx[same & idx$i != idx$j, ])] <- 1
  ins <- lv_insulation_score(m, window = 3)
  expect_length(ins, n)
  expect_true(is.na(ins[1]))
  expect_lt(ins[20], ins[10])
})

test_that("window contact counts sum pixels touching the boundary window", {
  bins <- data.frame(chrom = "chr1", start = seq(0, 9000, by = 1000),
                     end = seq(1000, 10000, by = 1000))
  pixels <- data.frame(bin1_id = c(1L, 4L, 9L), bin2_id = c(2L, 5L, 10L),
                       count = c(2, 3, 4))
  got <- lv_window_contacts(pixels$bin1_id, pixels$bin2_id, pixels$count,
                            bin_start = 4000, bin_size = 1000, flank = 1500)
  expect_equal(got, 3)
  # wider flank covers bins 1..10 entirely: all pixels counted
  got2 <- lv_window_contacts(pixels$bin1_id, pixels$bin2_id, pixels$count,
                             bin_start = 4000, bin_size = 1000, flank = 4500)
  expect_equal(got2, 9)
})
