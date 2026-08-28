test_that("read_cool parses bins and pixels", {
  tmp <- tempfile(fileext = ".cool")
  info <- make_test_cool(tmp)
  x <- read_cool(tmp)
  expect_s3_class(x, "lv_contacts")
  expect_equal(nrow(x$bins), nrow(info$bins))
  expect_equal(nrow(x$pixels), nrow(info$pixels))
  expect_true(all(x$pixels$count >= 1))
})

test_that("read_scool iterates cells", {
  tmp <- tempfile(fileext = ".scool")
  make_test_scool(tmp, n_cells = 3)
  cl <- read_scool(tmp, cells = c("cell1", "cell2"))
  expect_length(cl, 2)
  expect_equal(names(cl), c("cell1", "cell2"))
  expect_s3_class(cl[[1]], "lv_contacts")
})

test_that("read_mcool picks a resolution group", {
  tmp <- tempfile(fileext = ".mcool")
  skip_if_not(requireNamespace("hdf5r"))
  src <- tempfile(fileext = ".cool"); info <- make_test_cool(src)
  f <- hdf5r::H5File$new(tmp, mode = "w")
  g <- f$create_group("resolutions")$create_group(as.character(5000))
  n_bins <- nrow(info$bins)
  gc <- g$create_group("chroms"); gc[["name"]] <- "chr1"; gc[["length"]] <- n_bins * 5000L
  gb <- g$create_group("bins")
  gb[["chrom"]] <- info$bins$chrom; gb[["start"]] <- info$bins$start; gb[["end"]] <- info$bins$end
  gp <- g$create_group("pixels")
  gp[["bin1_id"]] <- info$pixels$bin1_id; gp[["bin2_id"]] <- info$pixels$bin2_id
  gp[["count"]] <- info$pixels$count
  f$close_all()
  x <- read_mcool(tmp, resolution = 5000)
  expect_s3_class(x, "lv_contacts")
  expect_equal(nrow(x$pixels), nrow(info$pixels))
})

test_that("read_pairs parses 4DN pairs incl. extra columns", {
  tmp <- tempfile()
  writeLines(c("## pairs format v1.0",
               "##columns: readID chr1 pos1 chr2 pos2 strand1 strand2",
               "r1\tchr1\t1000\tchr1\t8000\t+\t-",
               "r2\tchr1\t2000\tchr1\t9000\t+\t-\t1\t1"), tmp)
  x <- read_pairs(tmp)
  expect_equal(nrow(x), 2)
  expect_equal(x$pos1, c(1000L, 2000L))
  expect_equal(names(x)[8:9], c("extra1", "extra2"))
  expect_equal(x$extra2[2], 1)
})

test_that("read_sparse_contacts parses HiC-Pro triplets", {
  tmp <- tempfile()
  writeLines(c("0\t4000\t0\t4000\t3", "0\t4000\t8000\t12000\t1"), tmp)
  x <- read_sparse_contacts(tmp, chrom = "chr1", resolution = 4000)
  expect_equal(nrow(x$pixels), 2)
  expect_equal(x$pixels$count, c(3, 1))
})

test_that("import_boundaries reads BED", {
  tmp <- tempfile()
  writeLines(c("chr1\t50000\t55000\t20\tcellA\t1", "chr1\t50000\t55000\t25\tcellB\t1"), tmp)
  b <- import_boundaries(tmp)
  expect_equal(nrow(b), 2)
  expect_equal(b$contacts, c(20, 25))
})

test_that("read_hic fails cleanly without strawr", {
  skip_if_not(!requireNamespace("strawr", quietly = TRUE))
  expect_error(read_hic("x.hic", 5000, "chr1"), "strawr")
})
