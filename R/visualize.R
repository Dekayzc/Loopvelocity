# Velocity visualization on cell embeddings (supplementary notes,
# "extrapolation of cell state": flow from current to projected points).

# shared: Gaussian-weighted grid summary of per-cell arrows; returns and
# optionally draws grid arrows
.draw_grid_flow <- function(pos, arsd, grid.n, grid.sd, min.grid.cell.mass,
                            min.arrow.size, max.grid.arrow.length,
                            fixed.arrow.length, arrow.lwd, plot.grid.points) {
  rx <- range(c(range(pos[, 1]), range(pos[, 1] + arsd$xd)))
  ry <- range(c(range(pos[, 2]), range(pos[, 2] + arsd$yd)))
  gx <- seq(rx[1], rx[2], length.out = grid.n)
  gy <- seq(ry[1], ry[2], length.out = grid.n)
  if (is.null(grid.sd))
    grid.sd <- sqrt((gx[2] - gx[1])^2 + (gy[2] - gy[1])^2) / 2
  if (is.null(min.arrow.size))
    min.arrow.size <- sqrt((gx[2] - gx[1])^2 + (gy[2] - gy[1])^2) * 1e-2
  if (is.null(max.grid.arrow.length))
    max.grid.arrow.length <-
      sqrt(sum((graphics::par("pin") / c(length(gx), length(gy)))^2)) * 0.25

  garrows <- do.call(rbind, lapply(gx, function(x) {
    cd <- sqrt(outer(pos[, 2], -gy, "+")^2 + (x - pos[, 1])^2)
    cw <- stats::dnorm(cd, sd = grid.sd)
    gw <- Matrix::colSums(cw)
    cws <- pmax(1, Matrix::colSums(cw))
    gxd <- Matrix::colSums(cw * arsd$xd) / cws
    gyd <- Matrix::colSums(cw * arsd$yd) / cws
    al <- sqrt(gxd^2 + gyd^2)
    vg <- gw >= min.grid.cell.mass & al >= min.arrow.size
    cbind(rep(x, sum(vg)), gy[vg], x + gxd[vg], gy[vg] + gyd[vg])
  }))
  colnames(garrows) <- c("x0", "y0", "x1", "y1")
  if (nrow(garrows) > 0) {
    if (fixed.arrow.length) {
      suppressWarnings(graphics::arrows(garrows[, 1], garrows[, 2],
                                        garrows[, 3], garrows[, 4],
                                        length = 0.05, lwd = arrow.lwd))
    } else {
      pin <- graphics::par("pin"); usr <- graphics::par("usr")
      alen <- pmin(max.grid.arrow.length, sqrt(
        ((garrows[, 3] - garrows[, 1]) * pin[1] / diff(usr[c(1, 2)]))^2 +
        ((garrows[, 4] - garrows[, 2]) * pin[2] / diff(usr[c(3, 4)]))^2))
      suppressWarnings(apply(garrows, 1, function(x)
        graphics::arrows(x[1], x[2], x[3], x[4],
                         length = alen[1], lwd = arrow.lwd)))
    }
  }
  if (plot.grid.points)
    graphics::points(rep(gx, each = length(gy)), rep(gy, length(gx)),
                     pch = ".", cex = 1e-1, col = lv_alpha_color(1, alpha = 0.4))
  garrows
}

.cor_kernel <- function(emb, vel, n, corr.sigma, scale, n.cores, cell.colors,
                        cell.border.alpha, show.grid.flow, grid.n, grid.sd,
                        min.grid.cell.mass, min.arrow.size, arrow.scale,
                        max.grid.arrow.length, fixed.arrow.length,
                        plot.grid.points, arrow.lwd, xlab, ylab, do.par, ...) {
  if (do.par) graphics::par(mfrow = c(1, 1), mar = c(3.5, 3.5, 2.5, 1.5),
                            mgp = c(2, 0.65, 0), cex = 0.85)
  graphics::plot(emb, bg = cell.colors[rownames(emb)], pch = 21,
                 col = lv_alpha_color(1, alpha = cell.border.alpha),
                 xlab = xlab, ylab = ylab, ...)

  em <- as.matrix(vel$current)
  ccells <- intersect(rownames(emb), colnames(em))
  em <- em[, ccells, drop = FALSE]; emb <- emb[ccells, , drop = FALSE]
  nd <- as.matrix(vel$deltaE[, ccells, drop = FALSE])
  cgenes <- intersect(rownames(em), rownames(nd))
  nd <- nd[cgenes, , drop = FALSE]; em <- em[cgenes, , drop = FALSE]
  if (nrow(em) < 2)
    stop("velocity visualization needs at least 2 TADs in the velocity result")
  # keep only TADs present in deltaE (fitted rows); others would zero out
  # the delta correlations
  keep <- rownames(em) %in% rownames(vel$deltaE)
  if (sum(keep) < 2)
    stop("fewer than 2 TADs with identifiable steady-state fits")
  em <- em[keep, , drop = FALSE]; nd <- nd[keep, , drop = FALSE]

  cat("delta projections ... ")
  if (scale == "log") {
    cc <- lv_delta_cor_log10(em, (log10(abs(nd) + 1) * sign(nd)), nthreads = n.cores)
  } else if (scale == "sqrt") {
    cc <- lv_delta_cor_sqrt(em, (sqrt(abs(nd)) * sign(nd)), nthreads = n.cores)
  } else if (scale == "rank") {
    cc <- lv_delta_cor(apply(em, 2, rank), apply(nd, 2, rank), nthreads = n.cores)
  } else {
    cc <- lv_delta_cor(em, nd, nthreads = n.cores)
  }
  colnames(cc) <- rownames(cc) <- colnames(em)
  diag(cc) <- 0
  cat("done\n")

  cat("knn ... ")
  if (n > nrow(cc)) n <- nrow(cc)
  emb.knn <- balanced_knn(emb, k = n, maxl = nrow(emb), dist = "euclidean",
                          n.threads = n.cores)
  diag(emb.knn) <- 1
  cat("transition probs ... ")
  tp <- exp(cc / corr.sigma) * emb.knn
  tp <- t(t(tp) / Matrix::colSums(tp))
  tp <- as(tp, "dgCMatrix")
  cat("done\n")

  cat("calculating arrows ... ")
  arsd <- data.frame(t(lv_emb_arrows(emb, tp, arrow.scale, n.cores)))
  rownames(arsd) <- rownames(emb)
  ars <- data.frame(cbind(emb, emb + arsd))
  colnames(ars) <- c("x0", "y0", "x1", "y1")
  colnames(arsd) <- c("xd", "yd")
  cat("done\n")

  if (show.grid.flow) {
    cat("grid estimates ... ")
    garrows <- .draw_grid_flow(emb, arsd, grid.n, grid.sd, min.grid.cell.mass,
                               min.arrow.size, max.grid.arrow.length,
                               fixed.arrow.length, arrow.lwd, plot.grid.points)
    cat("done\n")
    return(invisible(list(tp = tp, cc = cc, garrows = garrows,
                          arrows = as.matrix(ars))))
  } else {
    pin <- graphics::par("pin"); usr <- graphics::par("usr")
    suppressWarnings(apply(ars, 1, function(x) {
      if (fixed.arrow.length) {
        graphics::arrows(x[1], x[2], x[3], x[4], length = 0.05, lwd = arrow.lwd)
      } else {
        ali <- sqrt(((x[3] - x[1]) * pin[1] / diff(usr[c(1, 2)]))^2 +
                    ((x[4] - x[2]) * pin[2] / diff(usr[c(3, 4)]))^2)
        graphics::arrows(x[1], x[2], x[3], x[4],
                         length = min(0.05, ali), lwd = arrow.lwd)
      }
    }))
    return(invisible(list(tp = tp, cc = cc, arrows = as.matrix(ars))))
  }
}

.eu_kernel <- function(emb, vel, n, sigma, beta, scale, embedding.knn,
                       control.for.neighborhood.density, n.cores, cell.colors,
                       cell.border.alpha, show.grid.flow, grid.n, grid.sd,
                       min.grid.cell.mass, min.arrow.size, arrow.scale,
                       max.grid.arrow.length, fixed.arrow.length,
                       plot.grid.points, arrow.lwd, xlab, ylab, do.par, ...) {
  if (do.par) graphics::par(mfrow = c(1, 1), mar = c(3.5, 3.5, 2.5, 1.5),
                            mgp = c(2, 0.65, 0), cex = 0.85)
  graphics::plot(emb, bg = cell.colors[rownames(emb)], pch = 21,
                 col = lv_alpha_color(1, alpha = cell.border.alpha),
                 xlab = xlab, ylab = ylab, ...)

  em <- as.matrix(vel$current); emn <- as.matrix(vel$projected)
  ccells <- intersect(rownames(emb), colnames(em))
  em <- em[, ccells, drop = FALSE]; emn <- emn[, ccells, drop = FALSE]; emb <- emb[ccells, , drop = FALSE]
  keep <- rownames(em) %in% rownames(vel$deltaE)
  em <- em[keep, , drop = FALSE]; emn <- emn[keep, , drop = FALSE]

  if (scale == "log") {
    em <- log10(em + 1); emn <- log10(emn + 1)
  } else if (scale == "sqrt") {
    em <- sqrt(em); emn <- sqrt(emn)
  }

  # cell-cell distances now and after extrapolation
  cc0 <- lv_col_euclid(em, em, nthreads = n.cores)
  cc <- lv_col_euclid(em, emn, nthreads = n.cores)
  diag(cc) <- diag(cc0)

  # automatic kernel scales from the distance distribution
  if (is.na(sigma)) sigma <- max(cc0[upper.tri(cc0)], na.rm = TRUE) / 10
  if (is.na(beta)) beta <- max(cc0[upper.tri(cc0)], na.rm = TRUE) / 20
  cat("sigma =", round(sigma, 3), " beta =", round(beta, 3),
      " transition probs ... ")
  f <- (cc / cc0)^beta; diag(f) <- 1
  tp <- exp(-((cc0 * f)^2) / (2 * sigma^2))
  np <- exp(-((cc0)^2) / (2 * sigma^2))
  colnames(tp) <- rownames(tp) <- colnames(np) <- rownames(np) <- colnames(em)

  if (n < nrow(emb)) {
    if (embedding.knn) {
      cell.knn <- balanced_knn(emb, k = n, maxl = nrow(emb), dist = "euclidean",
                               n.threads = n.cores)
    } else {
      cell.knn <- balanced_knn(t(em), k = n, maxl = ncol(em), dist = "cor",
                               n.threads = n.cores)
    }
    diag(cell.knn) <- 1
    tp <- tp * cell.knn; np <- np * cell.knn
  }
  tp <- t(t(tp) / Matrix::colSums(tp))
  np <- t(t(np) / Matrix::colSums(np))
  if (control.for.neighborhood.density) {
    np.f <- Matrix::diag(np)
    tp <- tp * np.f; np <- np * np.f
  }
  tp <- t(t(tp) / Matrix::colSums(tp))
  tp <- as(tp, "dgCMatrix")
  cat("done\n")

  cat("calculating arrows ... ")
  arsd <- data.frame(t(lv_emb_arrows(emb, tp, arrow.scale, n.cores)))
  rownames(arsd) <- rownames(emb)
  ars <- data.frame(cbind(emb, emb + arsd))
  colnames(ars) <- c("x0", "y0", "x1", "y1")
  colnames(arsd) <- c("xd", "yd")
  cat("done\n")

  if (show.grid.flow) {
    cat("grid estimates ... ")
    garrows <- .draw_grid_flow(emb, arsd, grid.n, grid.sd, min.grid.cell.mass,
                               min.arrow.size, max.grid.arrow.length,
                               fixed.arrow.length, arrow.lwd, plot.grid.points)
    cat("done\n")
    return(invisible(list(tp = tp, garrows = garrows, arrows = as.matrix(ars))))
  } else {
    pin <- graphics::par("pin"); usr <- graphics::par("usr")
    suppressWarnings(apply(ars, 1, function(x) {
      if (fixed.arrow.length) {
        graphics::arrows(x[1], x[2], x[3], x[4], length = 0.05, lwd = arrow.lwd)
      } else {
        ali <- sqrt(((x[3] - x[1]) * pin[1] / diff(usr[c(1, 2)]))^2 +
                    ((x[4] - x[2]) * pin[2] / diff(usr[c(3, 4)]))^2)
        graphics::arrows(x[1], x[2], x[3], x[4],
                         length = min(0.05, ali), lwd = arrow.lwd)
      }
    }))
    return(invisible(list(tp = tp, arrows = as.matrix(ars))))
  }
}

#' Plot loop velocity as a flow field on a cell embedding
#'
#' Draws the transition-weighted displacement of cells on a 2D embedding
#' (UMAP/tSNE/PCA): arrows point from the current chromatin contact state
#' toward the projected state. Two transition kernels are provided:
#' \code{"correlation"} (default; transition probabilities from
#' delta-correlations of the contact profiles within the embedding kNN) and
#' \code{"euclidean"} (transition probabilities from current-vs-projected
#' expression-space distances).
#'
#' @param emb matrix/data.frame of cells x 2 embedding coordinates with
#'   cell rownames
#' @param vel velocity result ([loop_velocity_estimates], optionally
#'   followed by [loop_velocity_em]/[project_cells])
#' @param kernel "correlation" or "euclidean" transition kernel
#' @param n neighbourhood size for the embedding kNN graph
#' @param cell.colors named colour vector for cells
#' @param corr.sigma kernel width for the correlation kernel
#' @param sigma,beta kernel parameters for the euclidean kernel (auto by
#'   default)
#' @param show.grid.flow summarize arrows on a Gaussian-weighted grid
#' @param grid.n grid points per axis
#' @param grid.sd Gaussian width of grid nodes
#' @param min.grid.cell.mass minimum summed kernel weight per node
#' @param min.arrow.size minimum plotted arrow length
#' @param max.grid.arrow.length arrow length cap when scaling
#' @param arrow.scale displacement scale multiplier
#' @param arrow.lwd arrow line width
#' @param fixed.arrow.length force uniform arrow heads
#' @param scale contact-space transform: "log" (default), "sqrt", "rank",
#'   "linear"
#' @param embedding.knn euclidean kernel: build the kNN in the embedding
#'   (TRUE) or in contact space (FALSE)
#' @param control.for.neighborhood.density euclidean kernel: compensate
#'   for cell-density variation
#' @param plot.grid.points mark grid nodes
#' @param n.cores worker threads for the C++ kernels
#' @param xlab,ylab axis labels
#' @param cell.border.alpha cell outline transparency
#' @param do.par reset graphics parameters before plotting
#' @param ... passed to \code{plot()}
#' @return invisibly, a list with the transition probability matrix
#'   \code{tp}, delta correlations \code{cc} (correlation kernel), arrow
#'   coordinates and grid arrows when requested
#' @examples
#' set.seed(1); n <- 30; th <- 0.7
#' # each TAD gets its own draw, with noise that gives every cell a
#' # residual from the constraint line -- the velocity signal
#' per_tad <- lapply(1:3, function(g) {
#'   Ig <- runif(n, 0.5, 3)
#'   list(I = Ig, S = (1 - th) / th * Ig + 2 + rnorm(n, 0, 0.05))
#' })
#' S <- t(sapply(per_tad, `[[`, "S")); I <- t(sapply(per_tad, `[[`, "I"))
#' dimnames(S) <- dimnames(I) <- list(paste0("tad", 1:3), paste0("c", 1:n))
#' vel <- loop_velocity_estimates(S, I, mode = "two-sided", kCells = 1,
#'                                mult = 1e9, verbose = FALSE)
#' emb <- cbind(runif(n), runif(n)); rownames(emb) <- colnames(S)
#' cc <- setNames(rep("grey60", n), colnames(S))
#' plot_loop_velocity(emb, vel, n = 10, cell.colors = cc,
#'                    show.grid.flow = TRUE, grid.n = 5)
#' @export
plot_loop_velocity <- function(emb, vel, kernel = c("correlation", "euclidean"),
                               n = 100, cell.colors = NULL, corr.sigma = 0.05,
                               sigma = NA, beta = 1, show.grid.flow = FALSE,
                               grid.n = 20, grid.sd = NULL,
                               min.grid.cell.mass = 1, min.arrow.size = NULL,
                               arrow.scale = 1, max.grid.arrow.length = NULL,
                               fixed.arrow.length = FALSE, plot.grid.points = FALSE,
                               scale = "log", embedding.knn = TRUE,
                               control.for.neighborhood.density = TRUE,
                               arrow.lwd = 1, xlab = "", ylab = "",
                               cell.border.alpha = 0.3,
                               n.cores = default_n_cores(), do.par = TRUE, ...) {
  kernel <- match.arg(kernel)
  if (is.null(cell.colors))
    cell.colors <- setNames(rep("white", nrow(emb)), rownames(emb))
  args <- list(emb = emb, vel = vel, n = n, cell.colors = cell.colors,
               show.grid.flow = show.grid.flow, grid.n = grid.n,
               grid.sd = grid.sd, min.grid.cell.mass = min.grid.cell.mass,
               min.arrow.size = min.arrow.size, arrow.scale = arrow.scale,
               max.grid.arrow.length = max.grid.arrow.length,
               fixed.arrow.length = fixed.arrow.length,
               plot.grid.points = plot.grid.points, arrow.lwd = arrow.lwd,
               xlab = xlab, ylab = ylab, do.par = do.par, ...)
  if (kernel == "correlation") {
    do.call(.cor_kernel, c(args, list(corr.sigma = corr.sigma, scale = scale,
                                      n.cores = n.cores,
                                      cell.border.alpha = cell.border.alpha)))
  } else {
    do.call(.eu_kernel, c(args, list(sigma = sigma, beta = beta, scale = scale,
                                     embedding.knn = embedding.knn,
                                     control.for.neighborhood.density =
                                       control.for.neighborhood.density,
                                     n.cores = n.cores,
                                     cell.border.alpha = cell.border.alpha)))
  }
}

#' Project velocity onto principal components
#'
#' Diagnostic plot: cells in consecutive PC pairs of the current contact
#' state, with arrows along the projected velocity increments.
#'
#' @param vel velocity result
#' @param nPcs number of successive PCs to visualize
#' @param cell.colors named colour vector
#' @param scale contact transform ("log" default)
#' @param plot.cols layout columns
#' @param norm.nPcs PCs for magnitude normalization (reserved)
#' @param arrow.scale displacement scale multiplier
#' @param arrow.lwd arrow line width
#' @param show.grid.flow,grid.n,grid.sd,min.grid.cell.mass,min.arrow.size
#'   grid-flow controls as in [plot_loop_velocity]
#' @param fixed.arrow.length,plot.grid.points,max.grid.arrow.length
#'   further grid rendering options
#' @param n.cores worker threads
#' @param do.par reset graphics parameters
#' @param ... passed to \code{plot()}
#' @param size.norm rescale current/projected by cell totals
#' @return invisibly, PCA info and projected velocity PCs
#' @examples
#' set.seed(1); n <- 25
#' th <- 0.7; Ig <- runif(n, 0.5, 3)
#' S <- matrix((1 - th) / th * Ig + 2, 3, n,
#'            dimnames = list(paste0("tad", 1:3), paste0("c", 1:n)))
#' I <- matrix(Ig, 3, n, dimnames = dimnames(S))
#' vel <- loop_velocity_estimates(S, I, mode = "two-sided", kCells = 1,
#'                                mult = 1e9, verbose = FALSE)
#' cc <- setNames(rep("grey60", n), colnames(S))
#' plot_velocity_pcs(vel, nPcs = 2, cell.colors = cc)
#' @export
plot_velocity_pcs <- function(vel, nPcs = 4, cell.colors = NULL, scale = "log",
                              plot.cols = min(3, nPcs - 1), norm.nPcs = NA,
                              show.grid.flow = FALSE, grid.n = 20,
                              grid.sd = NULL, min.grid.cell.mass = 1,
                              min.arrow.size = NULL, arrow.scale = 1,
                              arrow.lwd = 1, fixed.arrow.length = FALSE,
                              plot.grid.points = FALSE,
                              max.grid.arrow.length = NULL, size.norm = FALSE,
                              do.par = TRUE, n.cores = default_n_cores(), ...) {
  x <- as.matrix(vel$current); x1 <- as.matrix(vel$projected)
  if (scale == "log") { x <- log10(x + 1); x1 <- log10(x1 + 1) }
  else if (scale == "sqrt") { x <- sqrt(x); x1 <- sqrt(x1) }
  if (size.norm) {
    x <- t(t(x) / Matrix::colSums(x)); x1 <- t(t(x1) / Matrix::colSums(x1))
  }
  xc <- sweep(x, 1, Matrix::rowMeans(x))
  epc <- stats::prcomp(t(xc), center = FALSE, rank. = nPcs)
  nPcs <- min(nPcs, ncol(epc$x))
  scores <- epc$x                      # cells x PCs
  delta.pcs <- (t(x1) - t(x)) %*% epc$rotation   # cells x PCs
  delta.pcs <- delta.pcs * arrow.scale
  if (do.par) graphics::par(mfrow = c(ceiling((nPcs - 1) / plot.cols), plot.cols),
                            mar = c(3.5, 3.5, 2.5, 1.5), mgp = c(2, 0.65, 0),
                            cex = 0.85)
  invisible(lapply(seq_len(nPcs - 1), function(i) {
    pos <- scores[, c(i, i + 1), drop = FALSE]
    ppos <- pos + delta.pcs[, c(i, i + 1), drop = FALSE]
    graphics::plot(pos, bg = cell.colors[rownames(pos)], pch = 21,
                   col = lv_alpha_color(1, alpha = 0.3), lwd = 0.5,
                   xlab = paste("PC", i), ylab = paste("PC", i + 1),
                   main = paste("PC", i, "vs. PC", i + 1), ...)
    graphics::box()
    ars <- data.frame(pos[, 1], pos[, 2], ppos[, 1], ppos[, 2])
    colnames(ars) <- c("x0", "y0", "x1", "y1")
    arsd <- data.frame(xd = ars$x1 - ars$x0, yd = ars$y1 - ars$y0)
    rownames(ars) <- rownames(arsd) <- rownames(pos)
    if (show.grid.flow) {
      .draw_grid_flow(pos, arsd, grid.n, grid.sd, min.grid.cell.mass,
                      min.arrow.size, max.grid.arrow.length,
                      fixed.arrow.length, arrow.lwd, plot.grid.points)
    } else {
      suppressWarnings(apply(ars, 1, function(a)
        graphics::arrows(a[1], a[2], a[3], a[4], length = 0.05,
                         lwd = arrow.lwd)))
    }
  }))
  invisible(list(epc = epc, delta.pcs = delta.pcs))
}

#' tSNE view of the velocity field
#'
#' Fits a tSNE embedding on the (scaled) current contact state and draws
#' the velocity arrows on it via the correlation kernel. Requires the
#' \code{Rtsne} package.
#'
#' @param vel velocity result
#' @param cell.colors named colour vector
#' @param perplexity,nPcs passed to \code{Rtsne::Rtsne}
#' @param scale contact transform ("log" default)
#' @param ... passed to [plot_loop_velocity]
#' @return invisibly, the tSNE embedding and the plot details
#' @examples
#' \dontrun{
#' set.seed(1); n <- 50
#' th <- 0.7; Ig <- runif(n, 0.5, 3)
#' S <- matrix((1 - th) / th * Ig + 2, 3, n,
#'            dimnames = list(paste0("tad", 1:3), paste0("c", 1:n)))
#' I <- matrix(Ig, 3, n, dimnames = dimnames(S))
#' vel <- loop_velocity_estimates(S, I, mode = "two-sided", kCells = 1,
#'                                mult = 1e9, verbose = FALSE)
#' plot_velocity_tsne(vel, perplexity = 10)
#' }
#' @export
plot_velocity_tsne <- function(vel, cell.colors = NULL, perplexity = 30,
                               nPcs = 15, scale = "log", ...) {
  if (!requireNamespace("Rtsne", quietly = TRUE))
    stop("Rtsne is required for the tSNE view: install.packages('Rtsne')")
  x <- as.matrix(vel$current)
  if (scale == "log") x <- log10(x + 1)
  else if (scale == "sqrt") x <- sqrt(x)
  ts <- Rtsne::Rtsne(t(x), perplexity = perplexity, pca_center = TRUE,
                     initial_dims = nPcs, check_duplicates = FALSE)
  emb <- ts$Y; rownames(emb) <- colnames(x)
  details <- plot_loop_velocity(emb, vel, cell.colors = cell.colors, ...)
  invisible(list(emb = emb, details = details))
}

#' Keep TADs active in at least one cell cluster
#'
#' @param vel velocity result or list with element \code{S} (TAD x cell
#'   matrix)
#' @param clusters named vector mapping cells to cluster labels
#' @param min.max.cluster.average minimum cluster-mean contact level for a
#'   TAD to be kept
#' @return names of the retained TADs
#' @examples
#' m <- matrix(c(10, 0, 0, 0), 2, dimnames = list(c("t1", "t2"), c("a", "b")))
#' filter_tads_by_activity(list(S = m), clusters = c(a = "x", b = "y"))
#' @export
filter_tads_by_activity <- function(vel, clusters, min.max.cluster.average = 0.1) {
  S <- if (is.matrix(vel)) vel else vel$S
  clusters <- clusters[colnames(S)]
  cl <- sort(unique(clusters))
  cmeans <- sapply(cl, function(k) Matrix::rowMeans(S[, clusters == k, drop = FALSE]))
  if (is.null(dim(cmeans))) cmeans <- matrix(cmeans, nrow = nrow(S),
                                             dimnames = list(rownames(S), cl))
  keep <- apply(cmeans, 1, function(r) max(r) >= min.max.cluster.average)
  names(keep)[keep]
}

#' Diagnostics of a per-TAD steady-state fit
#'
#' Scatters size-normalized I- against S-boundary contacts for one TAD and
#' overlays the Eq 14 fit line (two-sided) or the Eq 13 slope (one-sided).
#'
#' @param vel output of [loop_velocity_estimates]
#' @param tad TAD name to show
#' @return invisibly, the fitted coefficients
#' @examples
#' set.seed(1); n <- 30
#' th <- 0.7; Ig <- runif(n, 0.5, 3)
#' S <- matrix((1 - th) / th * Ig + 2, 3, n,
#'            dimnames = list(paste0("tad", 1:3), paste0("c", 1:n)))
#' I <- matrix(Ig, 3, n, dimnames = dimnames(S))
#' vel <- loop_velocity_estimates(S, I, mode = "two-sided", kCells = 1,
#'                                mult = 1e9, verbose = FALSE)
#' plot_si_fit(vel, "tad1")
#' @export
plot_si_fit <- function(vel, tad) {
  s <- as.numeric(vel$S[tad, ]); i <- as.numeric(vel$I[tad, ])
  graphics::plot(i, s, pch = 21, bg = lv_alpha_color(1, 0.3),
                 xlab = "I boundary contacts", ylab = "S boundary contacts",
                 main = tad)
  if (vel$mode[tad] == "two-sided") {
    a <- (1 - vel$theta[tad]) / vel$theta[tad]
    b <- if (is.na(vel$gamma[tad])) 0 else (1 / vel$gamma[tad]) * (2 - 1 / vel$theta[tad])
    graphics::abline(a = b, b = a, lty = 2, col = 2)
    graphics::legend("bottomright", bty = "n",
                     legend = sprintf("theta = %.3f\ngamma = %.3f",
                                      vel$theta[tad], vel$gamma[tad]))
    invisible(c(intercept = unname(b), slope = unname(a)))
  } else {
    graphics::abline(a = 0, b = vel$ratio[tad], lty = 2, col = 2)
    graphics::legend("bottomright", bty = "n",
                     legend = sprintf("gamma/beta = %.3f", vel$ratio[tad]))
    invisible(c(intercept = 0, slope = unname(vel$ratio[tad])))
  }
}
