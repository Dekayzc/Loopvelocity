# Hi-C contact readers: cooler (.cool/.mcool/.scool), 4DN pairs, sparse
# triplet matrices, and boundary BED tables.

new_lv_contacts <- function(bins, pixels) {
  structure(list(bins = bins, pixels = pixels), class = "lv_contacts")
}

.cool_group_read <- function(h5group) {
  bins <- data.frame(
    chrom = as.character(h5group[["bins/chrom"]][]),
    start = as.numeric(h5group[["bins/start"]][]),
    end = as.numeric(h5group[["bins/end"]][]))
  pixels <- data.frame(
    bin1_id = as.integer(h5group[["pixels/bin1_id"]][]) + 1L,
    bin2_id = as.integer(h5group[["pixels/bin2_id"]][]) + 1L,
    count = as.numeric(h5group[["pixels/count"]][]))
  new_lv_contacts(bins, pixels)
}

.require_hdf5r <- function() {
  if (!requireNamespace("hdf5r", quietly = TRUE))
    stop("hdf5r is required for cooler-format input: install.packages('hdf5r') ",
         "or use the conda-forge r-hdf5r package")
}

#' Read a single-resolution cooler contact file
#'
#' Parses the cooler layout (\code{chroms}, \code{bins}, \code{pixels}
#' groups) into a per-cell contact record.
#'
#' @param file path to a .cool file
#' @return object of class \code{lv_contacts}: a list with \code{bins}
#'   (data.frame: chrom, start, end) and \code{pixels} (data.frame: bin1_id,
#'   bin2_id, count; 1-based bin indices)
#' @examples
#' if (requireNamespace("hdf5r")) {
#'   f <- tempfile(fileext = ".cool")
#'   h <- hdf5r::H5File$new(f, mode = "w")
#'   g <- h$create_group("chroms"); g[["name"]] <- "chr1"; g[["length"]] <- 5e4
#'   g <- h$create_group("bins")
#'   g[["chrom"]] <- rep("chr1", 5)
#'   g[["start"]] <- seq(0, 4e4, by = 1e4); g[["end"]] <- seq(1e4, 5e4, by = 1e4)
#'   g <- h$create_group("pixels")
#'   g[["bin1_id"]] <- 0:2; g[["bin2_id"]] <- 1:3; g[["count"]] <- c(2, 1, 3)
#'   h$close_all()
#'   x <- read_cool(f)
#'   head(x$pixels)
#' }
#' @export
read_cool <- function(file) {
  .require_hdf5r()
  f <- hdf5r::H5File$new(file, mode = "r"); on.exit(f$close_all(), add = TRUE)
  .cool_group_read(f)
}

#' Read one resolution layer of a multi-resolution cooler file
#'
#' @param file path to a .mcool file
#' @param resolution resolution (bp) selecting the group under
#'   \code{/resolutions}
#' @return \code{lv_contacts}
#' @examples
#' if (requireNamespace("hdf5r")) {
#'   f <- tempfile(fileext = ".mcool")
#'   h <- hdf5r::H5File$new(f, mode = "w")
#'   g <- h$create_group("resolutions")$create_group("10000")
#'   gc <- g$create_group("chroms"); gc[["name"]] <- "chr1"; gc[["length"]] <- 5e4
#'   gb <- g$create_group("bins")
#'   gb[["chrom"]] <- rep("chr1", 5)
#'   gb[["start"]] <- seq(0, 4e4, by = 1e4); gb[["end"]] <- seq(1e4, 5e4, by = 1e4)
#'   gp <- g$create_group("pixels")
#'   gp[["bin1_id"]] <- 0:1; gp[["bin2_id"]] <- 1:2; gp[["count"]] <- c(2, 1)
#'   h$close_all()
#'   x <- read_mcool(f, resolution = 10000)
#'   head(x$pixels)
#' }
#' @export
read_mcool <- function(file, resolution) {
  .require_hdf5r()
  f <- hdf5r::H5File$new(file, mode = "r"); on.exit(f$close_all(), add = TRUE)
  g <- f[[paste0("resolutions/", as.character(resolution))]]
  .cool_group_read(g)
}

#' Read a single-cell cooler collection (.scool)
#'
#' Iterates \code{/cells/<cell>} groups, each holding a complete cooler
#' layout for one cell.
#'
#' @param file path to a .scool file
#' @param cells optional subset of cell names; defaults to all
#' @return named list of \code{lv_contacts}, one per cell
#' @examples
#' if (requireNamespace("hdf5r")) {
#'   f <- tempfile(fileext = ".scool")
#'   h <- hdf5r::H5File$new(f, mode = "w")
#'   cg <- h$create_group("cells")
#'   for (cn in c("c1", "c2")) {
#'     g <- cg$create_group(cn)
#'     gc <- g$create_group("chroms"); gc[["name"]] <- "chr1"; gc[["length"]] <- 3e4
#'     gb <- g$create_group("bins")
#'     gb[["chrom"]] <- rep("chr1", 3)
#'     gb[["start"]] <- c(0, 1e4, 2e4); gb[["end"]] <- c(1e4, 2e4, 3e4)
#'     gp <- g$create_group("pixels")
#'     gp[["bin1_id"]] <- 0:1; gp[["bin2_id"]] <- 1:2; gp[["count"]] <- c(1, 2)
#'   }
#'   h$close_all()
#'   cl <- read_scool(f, cells = "c1")
#'   cl$c1$pixels
#' }
#' @export
read_scool <- function(file, cells = NULL) {
  .require_hdf5r()
  f <- hdf5r::H5File$new(file, mode = "r"); on.exit(f$close_all(), add = TRUE)
  avail <- hdf5r::list.groups(f[["cells"]], recursive = FALSE)
  if (is.null(cells)) cells <- avail
  if (!all(cells %in% avail))
    stop("unknown cells: ", paste(setdiff(cells, avail), collapse = ", "))
  stats::setNames(lapply(cells, function(cn)
    .cool_group_read(f[[paste0("cells/", cn)]])), cells)
}

#' Read a 4DN pairs file
#'
#' Tolerant to extra columns beyond the seven standard ones (mapq and
#' similar annotations are kept as \code{extra1}, \code{extra2}, ...).
#'
#' @param file path to a .pairs or .pairs.gz file (header lines starting
#'   with \code{#} are skipped)
#' @return data.frame with columns readID, chr1, pos1, chr2, pos2,
#'   strand1, strand2, plus any extra columns
#' @examples
#' f <- tempfile()
#' writeLines(c("## pairs format v1.0",
#'   "r1\\tchr1\\t1000\\tchr1\\t8000\\t+\\t-\\t1\\t1"), f)
#' p <- read_pairs(f)
#' p$pos1
#' @export
read_pairs <- function(file) {
  ln <- utils::read.delim(file, header = FALSE, comment.char = "#")
  std <- c("readID", "chr1", "pos1", "chr2", "pos2", "strand1", "strand2")
  if (ncol(ln) < length(std))
    stop("pairs file has ", ncol(ln), " columns; expected at least 7")
  names(ln) <- c(std, paste0("extra", seq_len(ncol(ln) - length(std))))
  ln
}

#' Read a sparse contact-triplet matrix (HiC-Pro / Homer style)
#'
#' @param file path to the matrix file: rows are
#'   \code{start1 end1 start2 end2 count}
#' @param chrom chromosome name to attach to the bins
#' @param resolution bin size (bp) used to discretize positions
#' @return \code{lv_contacts}
#' @examples
#' f <- tempfile()
#' writeLines(c("0 10000 0 10000 3", "0 10000 20000 30000 1"), f)
#' x <- read_sparse_contacts(f, chrom = "chr1", resolution = 10000)
#' x$pixels
#' @export
read_sparse_contacts <- function(file, chrom, resolution) {
  # sep = "": accepts both tab- and space-separated triplet files
  ln <- utils::read.delim(file, header = FALSE, sep = "",
                          col.names = c("start1", "end1", "start2", "end2", "count"))
  pix <- data.frame(
    bin1_id = floor(ln$start1 / resolution) + 1L,
    bin2_id = floor(ln$start2 / resolution) + 1L,
    count = ln$count)
  nmax <- max(pix$bin1_id, pix$bin2_id)
  bins <- data.frame(chrom = chrom,
                     start = (seq_len(nmax) - 1) * resolution,
                     end = seq_len(nmax) * resolution)
  new_lv_contacts(bins, pix)
}

#' Import a TAD-boundary contact table
#'
#' Reads the boundary summary table produced by upstream TAD callers:
#' columns chr, start, end, contacts, cell, regionid. This is the format
#' consumed by [read_tad_contacts].
#'
#' @param bed path to the boundary table
#' @return data.frame
#' @examples
#' f <- tempfile()
#' writeLines(c("chr1\t50000\t55000\t20\tcellA\t1",
#'             "chr1\t50000\t55000\t25\tcellB\t1"), f)
#' import_boundaries(f)
#' @export
import_boundaries <- function(bed) {
  utils::read.delim(bed, header = FALSE,
                    col.names = c("chr", "start", "end", "contacts", "cell", "regionid"),
                    stringsAsFactors = FALSE)
}
