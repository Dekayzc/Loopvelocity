# Juicer .hic contact extraction (optional, via strawr).

#' Read contacts from a Juicer .hic file
#'
#' @param file path to a .hic file
#' @param resolution resolution (bp) at which to pull contacts
#' @param chroms chromosomes to extract, e.g. \code{c("chr1", "chr2")};
#'   cis-contacts per chromosome
#' @return data.frame of contact pixels (chr1, x1, chr2, y1, count)
#' @examples
#' \dontrun{
#' px <- read_hic("allpairs.hic", resolution = 5000, chroms = c("chr1", "chr2"))
#' }
#' @export
read_hic <- function(file, resolution, chroms = NULL) {
  if (!requireNamespace("strawr", quietly = TRUE))
    stop("strawr is required for Juicer .hic input: install.packages('strawr')")
  if (is.null(chroms)) stop("specify chroms, e.g. c('chr1','chr2')")
  do.call(rbind, lapply(chroms, function(ch)
    as.data.frame(strawr::straw(paste0("file://", file), ch, ch, "BP", resolution))))
}
