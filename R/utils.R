# Shared utilities: cell pooling graphs, colour helpers, defaults.

#' Balanced k-nearest-neighbour graph on cells
#'
#' Builds a balanced kNN graph over single cells from their contact profiles
#' and normalizes edge weights to a fixed total per cell. Neighbourhood
#' balancing prevents highly similar cells from dominating the pooled
#' estimates. Loop velocity follows the velocity pooling framework for this
#' step (supplementary notes, "Inference of parameters": balanced KNN
#' pooling).
#'
#' @param x cells x features matrix (e.g. log-normalized boundary contacts)
#' @param k number of neighbours
#' @param maxl cap on how often a cell may be chosen as a neighbour
#' @param return.distance.values keep edge distances instead of unit weights
#' @param n.threads OpenMP threads for the sort phase
#' @param dist distance flavour: "cor" (1 - Pearson) or "euclidean";
#'   alternatively pass a \code{dist} object directly
#' @return dense cell x cell adjacency matrix with rownames/colnames from x
#' @examples
#' set.seed(1)
#' m <- matrix(rnorm(60), nrow = 20, dimnames = list(paste0("c", 1:20), NULL))
#' knn <- balanced_knn(m, k = 5)
#' dim(knn); rowSums(knn)
#' @export
balanced_knn <- function(x, k, maxl = k, return.distance.values = FALSE,
                         n.threads = 1, dist = "cor") {
  if (inherits(dist, "dist")) {
    if (!all(labels(dist) == rownames(x)))
      stop("balanced_knn(): supplied distance does not match the rows of x")
    cd <- as.matrix(dist)
  } else if (dist == "cor") {
    cd <- 1 - lv_col_cor(t(x))
  } else if (dist == "euclidean") {
    cd <- as.matrix(stats::dist(x))
  } else {
    stop(paste("unknown distance", dist, "specified"))
  }
  z <- lv_balanced_knn_graph(cd, k, maxl, return.distance.values, n.threads)
  z <- as.matrix(z)
  rownames(z) <- colnames(z) <- rownames(x)
  z
}

#' Set a vector's names to its own values
#'
#' @param x vector
#' @return x with \code{names(x) <- x}
#' @examples
#' lv_set_names(c("a", "b"))
#' @export
lv_set_names <- function(x) { stats::setNames(x, x) }

#' Adjust colour transparency, keeping vector names
#'
#' @param x colour vector
#' @param alpha transparency passed to \code{adjustcolor} as \code{alpha.f}
#' @param ... further arguments to \code{adjustcolor}
#' @return colour vector with original names
#' @examples
#' lv_alpha_color(c("red", "blue"), alpha = 0.5)
#' @export
lv_alpha_color <- function(x, alpha = 1, ...) {
  y <- grDevices::adjustcolor(x, alpha.f = alpha, ...)
  names(y) <- names(x)
  y
}

#' Map numeric values onto a colour gradient
#'
#' @param x numeric vector; nonnegative values use a sequential ramp,
#'   signed values a diverging blue-grey-red ramp
#' @param gradientPalette optional explicit gradient
#' @param zlim optional value range
#' @param gradient.range.quantile quantile clip for automatic zlim
#' @return colour vector (named if x was)
#' @examples
#' lv_val2col(c(-2, -1, 0, 1, 2))
#' @export
lv_val2col <- function(x, gradientPalette = NULL, zlim = NULL,
                       gradient.range.quantile = 0.95) {
  if (all(sign(x) >= 0)) {
    if (is.null(gradientPalette))
      gradientPalette <- grDevices::colorRampPalette(c("gray90", "red"), space = "Lab")(1024)
    if (is.null(zlim)) {
      zlim <- as.numeric(stats::quantile(na.omit(x),
        p = c(1 - gradient.range.quantile, gradient.range.quantile)))
      if (diff(zlim) == 0) zlim <- as.numeric(range(na.omit(x)))
    }
  } else {
    if (is.null(gradientPalette))
      gradientPalette <- grDevices::colorRampPalette(c("blue", "grey90", "red"), space = "Lab")(1024)
    if (is.null(zlim)) {
      zlim <- c(-1, 1) * as.numeric(stats::quantile(na.omit(abs(x)), p = gradient.range.quantile))
      if (diff(zlim) == 0) zlim <- c(-1, 1) * as.numeric(na.omit(max(abs(x))))
    }
  }
  x[x < zlim[1]] <- zlim[1]; x[x > zlim[2]] <- zlim[2]
  x <- (x - zlim[1]) / (zlim[2] - zlim[1])
  gp <- gradientPalette[x * (length(gradientPalette) - 1) + 1]
  if (!is.null(names(x))) names(gp) <- names(x)
  gp
}

#' Default contact-intensity palette for embedding overlays
#' @return vector of 256 colours
#' @examples
#' plot(1:10, col = lv_expression_palette()[seq(1, 256, length.out = 10)], pch = 19)
#' @export
lv_expression_palette <- function() {
  grDevices::colorRampPalette(
    c("grey90", "yellow2", "orange2", "red3", "darkred"))(256)
}

#' Default worker count (physical cores)
#' @examples
#' default_n_cores()
#' @export
default_n_cores <- function() max(1L, parallel::detectCores(logical = FALSE))
