# Test fixtures: fabricate minimal cooler-layout HDF5 files.

make_test_cool <- function(path, chrom = "chr1", n_bins = 10, n_pixels = 30, seed = 1) {
  if (!requireNamespace("hdf5r")) skip("hdf5r not available")
  set.seed(seed)
  bins <- data.frame(chrom = chrom,
                     start = seq(0, by = 5000, length.out = n_bins),
                     end = seq(5000, by = 5000, length.out = n_bins))
  i <- sample(n_bins, n_pixels, replace = TRUE)
  j <- pmin(n_bins, i + sample(0:3, n_pixels, replace = TRUE))
  keep <- i < j
  pixels <- data.frame(bin1_id = i[keep] - 1L, bin2_id = j[keep] - 1L,
                       count = sample(1:5, sum(keep), replace = TRUE))
  f <- hdf5r::H5File$new(path, mode = "w")
  on.exit(f$close_all(), add = TRUE)
  g <- f$create_group("chroms"); g[["name"]] <- chrom; g[["length"]] <- n_bins * 5000L
  g <- f$create_group("bins")
  g[["chrom"]] <- bins$chrom; g[["start"]] <- bins$start; g[["end"]] <- bins$end
  g <- f$create_group("pixels")
  g[["bin1_id"]] <- pixels$bin1_id; g[["bin2_id"]] <- pixels$bin2_id; g[["count"]] <- pixels$count
  invisible(list(bins = bins, pixels = pixels))
}

make_test_scool <- function(path, n_cells = 3) {
  if (!requireNamespace("hdf5r")) skip("hdf5r not available")
  set.seed(0)
  n_bins <- 10
  bins <- data.frame(chrom = "chr1", start = seq(0, by = 5000, length.out = n_bins),
                     end = seq(5000, by = 5000, length.out = n_bins))
  f <- hdf5r::H5File$new(path, mode = "w")
  on.exit(f$close_all(), add = TRUE)
  cg <- f$create_group("cells")
  for (cn in paste0("cell", seq_len(n_cells))) {
    g <- cg$create_group(cn)
    gc <- g$create_group("chroms"); gc[["name"]] <- "chr1"; gc[["length"]] <- n_bins * 5000L
    gb <- g$create_group("bins")
    gb[["chrom"]] <- bins$chrom; gb[["start"]] <- bins$start; gb[["end"]] <- bins$end
    gp <- g$create_group("pixels")
    i <- sample(n_bins - 1, 20, TRUE)
    j <- pmin(n_bins, i + sample(1:3, 20, TRUE))
    gp[["bin1_id"]] <- i - 1L; gp[["bin2_id"]] <- j - 1L; gp[["count"]] <- sample(1:5, 20, TRUE)
  }
  invisible(f)
}
