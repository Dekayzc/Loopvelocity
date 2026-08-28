test_that("call_tad_boundaries finds the engineered boundary", {
  n <- 40; m <- matrix(0, n, n)
  idx <- expand.grid(i = 1:n, j = 1:n)
  same <- (idx$i <= 20) == (idx$j <= 20)
  m[as.matrix(idx[same & idx$i != idx$j, ])] <- 1
  b <- call_tad_boundaries(m, window = 3, min_delta = 0.1)
  expect_equal(nrow(b), 1)
  expect_equal(b$bin, 20)
})

test_that("assign_si_boundaries swaps high-contact side to S and tests", {
  set.seed(4); n_tads <- 30; n_cells <- 40
  L <- matrix(rpois(n_tads * n_cells, 8), n_tads)
  R <- matrix(rpois(n_tads * n_cells, 3), n_tads)
  rownames(L) <- rownames(R) <- paste0("tad", seq_len(n_tads))
  colnames(L) <- colnames(R) <- paste0("c", seq_len(n_cells))
  res <- assign_si_boundaries(L, R, test = "t", padj = 0.05)
  expect_true(all(rowSums(res$S) >= rowSums(res$I)))
  expect_true(all(res$mode %in% c("one-sided", "two-sided")))
  # clearly different contact levels -> significant -> two-sided
  expect_equal(as.character(res$mode), rep("two-sided", n_tads))
})

test_that("read_tad_contacts builds S/I matrices from BEDs", {
  tmp_l <- tempfile(); tmp_r <- tempfile()
  # left/right boundaries of a TAD sit at different coordinates; the
  # regionid enumerates the TAD in both tables
  writeLines(c("chr1\t50000\t55000\t10\tcellA\t1", "chr1\t50000\t55000\t12\tcellB\t1",
               "chr1\t90000\t95000\t3\tcellA\t2", "chr1\t90000\t95000\t1\tcellB\t2"), tmp_l)
  writeLines(c("chr1\t120000\t125000\t2\tcellA\t1", "chr1\t120000\t125000\t1\tcellB\t1",
               "chr1\t160000\t165000\t8\tcellA\t2", "chr1\t160000\t165000\t9\tcellB\t2"), tmp_r)
  out <- read_tad_contacts(tmp_l, tmp_r, assign = FALSE)
  expect_equal(dim(out$left), c(2, 2))
  expect_equal(dim(out$right), c(2, 2))
  expect_equal(rownames(out$left), rownames(out$right))
  expect_equal(out$left["chr1-1", "cellA"], 10)
  expect_equal(out$right["chr1-2", "cellB"], 9)
})
