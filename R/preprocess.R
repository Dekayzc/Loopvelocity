# From raw contacts to S/I boundary matrices: insulation-based boundary
# calling, per-cell boundary-window contact extraction, S/I assignment and
# one/two-sided extrusion mode testing (supplementary notes, "Determining
# whether the loop extrusion process is one-sided or two-sided").

#' Diamond insulation score of a contact matrix
#'
#' Thin wrapper over the C++ implementation; see \code{\link{lv_insulation_score}}.
#'
#' @param m square dense contact matrix
#' @param window diamond radius in bins
#' @examples
#' # two contact domains meeting at bin 5: insulation dips at the boundary
#' n <- 20; m <- matrix(0, n, n)
#' idx <- expand.grid(i = 1:n, j = 1:n)
#' same <- (idx$i <= 5) == (idx$j <= 5)
#' m[as.matrix(idx[same & idx$i != idx$j, ])] <- 1
#' ins <- compute_insulation_score(m, window = 3)
#' plot(ins, type = "l"); abline(v = 5, lty = 2)
#' @export
compute_insulation_score <- function(m, window) {
  stopifnot(is.matrix(m), nrow(m) == ncol(m))
  lv_insulation_score(m, window = window)
}

#' Call TAD boundaries from an insulation profile
#'
#' A bin is a boundary when its insulation score is a local minimum and the
#' dip relative to the flanking minima (the delta) falls below the
#' \code{min_delta} quantile of the delta distribution.
#'
#' @param m square dense contact matrix
#' @param window diamond radius for insulation
#' @param min_delta lower-quantile cutoff on the insulation dip
#' @return data.frame with columns bin, insulation, delta
#' @examples
#' n <- 40; m <- matrix(0, n, n)
#' idx <- expand.grid(i = 1:n, j = 1:n)
#' same <- (idx$i <= 20) == (idx$j <= 20)
#' m[as.matrix(idx[same & idx$i != idx$j, ])] <- 1
#' call_tad_boundaries(m, window = 3, min_delta = 0.1)
#' @export
call_tad_boundaries <- function(m, window, min_delta = 0.2) {
  ins <- compute_insulation_score(m, window)
  n <- length(ins)
  delta <- rep(NA_real_, n)
  local_min <- rep(FALSE, n)
  for (i in (window + 1):(n - window)) {
    left <- ins[max(1, i - window):(i - 1)]
    right <- ins[(i + 1):min(n, i + window)]
    if (any(is.na(left)) || any(is.na(right))) next
    delta[i] <- ins[i] - (min(left) + min(right)) / 2
    local_min[i] <- ins[i] < ins[i - 1] && ins[i] <= ins[i + 1]
  }
  thr <- stats::quantile(delta, probs = min_delta, na.rm = TRUE)
  keep <- which(!is.na(delta) & local_min & delta <= min(thr, 0))
  data.frame(bin = keep, insulation = ins[keep], delta = delta[keep])
}

#' Per-cell contacts around TAD boundaries
#'
#' For each cell and each boundary, sums contact pixels with at least one
#' endpoint inside the boundary window [start - flank, end + flank).
#'
#' @param cell_contacts named list of \code{lv_contacts} (one per cell, same
#'   bin grid)
#' @param boundaries data.frame with columns chrom, start, end on that grid
#' @param flank window extension in bp (default 25 kb -> ~50 kb windows)
#' @return TAD x cell matrix of window contact counts
#' @examples
#' bins <- data.frame(chrom = "chr1", start = c(0, 1e4, 2e4), end = c(1e4, 2e4, 3e4))
#' px <- data.frame(bin1_id = c(1L, 1L, 2L), bin2_id = c(2L, 3L, 3L), count = c(2, 1, 4))
#' cells <- list(cA = structure(list(bins = bins, pixels = px), class = "lv_contacts"))
#' bnd <- data.frame(chrom = "chr1", start = 1e4, end = 2e4)
#' boundary_contacts(cells, bnd, flank = 5000)
#' @export
boundary_contacts <- function(cell_contacts, boundaries, flank = 25e3) {
  stopifnot(is.list(cell_contacts), nrow(boundaries) >= 1)
  ref_bins <- cell_contacts[[1]]$bins
  bin_size <- ref_bins$end[1] - ref_bins$start[1]
  res <- matrix(0, nrow = nrow(boundaries), ncol = length(cell_contacts))
  rownames(res) <- paste(boundaries$chrom, boundaries$start,
                         seq_len(nrow(boundaries)), sep = "-")
  colnames(res) <- names(cell_contacts)
  for (ci in seq_along(cell_contacts)) {
    px <- cell_contacts[[ci]]$pixels
    bn <- cell_contacts[[ci]]$bins
    # map each boundary to the bin index on this cell's grid
    bin_idx <- vapply(seq_len(nrow(boundaries)), function(bi)
      which(bn$chrom == boundaries$chrom[bi] &
            bn$start <= boundaries$start[bi] & bn$end > boundaries$start[bi]),
      integer(1))
    for (bi in seq_len(nrow(boundaries))) {
      if (length(bin_idx[bi]) == 0 || is.na(bin_idx[bi])) next
      res[bi, ci] <- lv_window_contacts(px$bin1_id, px$bin2_id, px$count,
                                        bin_start = bn$start[bin_idx[bi]],
                                        bin_size = bin_size, flank = flank)
    }
  }
  res
}

#' Assign S/I boundary roles and extrusion mode per TAD
#'
#' The boundary with more contacts in a TAD is the S (stabilizing)
#' boundary, the other the I (initiation) boundary. A per-TAD significance
#' test of the two boundaries' contact numbers decides the extrusion mode:
#' no significant difference -> one-sided, significant -> two-sided
#' (supplementary notes, mode determination).
#'
#' @param left,right TAD x cell matrices of contacts around the two
#'   boundaries of each TAD
#' @param test "t" (default, per the notes) or "wilcox" (rank-based)
#' @param padj adjusted-p cutoff for the two-sided call (BH correction)
#' @return list of class \code{lv_si}: \code{S}, \code{I} matrices,
#'   \code{mode} vector, \code{p.value}, \code{padj}, \code{swapped}
#' @examples
#' set.seed(1)
#' L <- matrix(rpois(30 * 40, 8), 30, dimnames = list(paste0("tad", 1:30), paste0("c", 1:40)))
#' R <- matrix(rpois(30 * 40, 3), 30, dimnames = dimnames(L))
#' si <- assign_si_boundaries(L, R, test = "t", padj = 0.05)
#' table(si$mode)          # clearly different levels -> two-sided
#' @export
assign_si_boundaries <- function(left, right, test = c("t", "wilcox"), padj = 0.05) {
  test <- match.arg(test)
  stopifnot(all(dim(left) == dim(right)), all(colnames(left) == colnames(right)))
  n <- nrow(left)
  S <- left; I <- right
  swap <- as.vector(Matrix::rowSums(right) > Matrix::rowSums(left))
  if (any(swap)) {
    S[swap, ] <- right[swap, , drop = FALSE]
    I[swap, ] <- left[swap, , drop = FALSE]
  }
  p <- vapply(seq_len(n), function(i)
    tryCatch(
      if (test == "t") stats::t.test(left[i, ], right[i, ])$p.value
      else stats::wilcox.test(left[i, ], right[i, ])$p.value,
      error = function(e) NA_real_),
    numeric(1))
  padj_v <- stats::p.adjust(p, method = "BH")
  out <- list(S = S, I = I,
              mode = ifelse(!is.na(padj_v) & padj_v < padj,
                            "two-sided", "one-sided"),
              p.value = p, padj = padj_v, swapped = swap)
  structure(out, class = "lv_si")
}

#' Read boundary contact tables into S/I matrices
#'
#' Quick path for pre-summarized boundary tables (chr, start, end,
#' contacts, cell, regionid): pivots them to TAD x cell matrices and runs
#' [assign_si_boundaries].
#'
#' @param left_bed,right_bed paths to the two boundary tables
#' @param assign run S/I assignment (default) or return raw left/right
#' @param test,padj passed to [assign_si_boundaries]
#' @return \code{lv_si} list, or \code{list(left, right)} when assign = FALSE
#' @examples
#' fl <- tempfile(); fr <- tempfile()
#' writeLines(c("chr1\t50000\t55000\t10\tcellA\t1",
#'             "chr1\t50000\t55000\t12\tcellB\t1"), fl)
#' writeLines(c("chr1\t120000\t125000\t2\tcellA\t1",
#'             "chr1\t120000\t125000\t1\tcellB\t1"), fr)
#' si <- read_tad_contacts(fl, fr)   # left/right pair by chr + region id
#' si$S[, c("cellA", "cellB")]
#' @export
read_tad_contacts <- function(left_bed, right_bed, assign = TRUE,
                              test = c("t", "wilcox"), padj = 0.05) {
  lb <- import_boundaries(left_bed); rb <- import_boundaries(right_bed)
  .long_to_mat <- function(df) {
    # regionid enumerates the TADs of a chromosome per cell and is shared
    # between the left/right boundary tables; the boundary START differs
    # between the two sides of a TAD, so chr+regionid (not start) is the
    # pairing key
    df$key <- paste(df$chr, df$regionid, sep = "-")
    cells <- sort(unique(df$cell)); regions <- sort(unique(df$key))
    m <- matrix(0, length(regions), length(cells), dimnames = list(regions, cells))
    m[cbind(match(df$key, regions), match(df$cell, cells))] <- df$contacts
    m
  }
  left <- .long_to_mat(lb); right <- .long_to_mat(rb)
  if (assign) assign_si_boundaries(left, right, test = test, padj = padj)
  else list(left = left, right = right)
}
