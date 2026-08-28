# Steady-state loop velocity estimation (supplementary notes,
# "Steady-state model"). Two-sided: per TAD, fit S ~ a*I + b (Eq 14 with
# alpha == 1) and solve theta = 1/(1+a), gamma = (2 - 1/theta)/b.
# One-sided: fit S ~ a*I + b (Eq 13) and report ratio = 1/a = gamma/beta
# (scale-free steady-state constraint).

#' Steady-state loop velocity estimates
#'
#' Fits the steady-state constraint of the loop extrusion model to per-TAD
#' S/I boundary contact matrices. Balanced kNN pooling over cells and
#' top/bottom quantile weighting follow the velocity framework referenced
#' in the supplementary notes.
#'
#' Two-sided parameterization phi(alpha, theta, gamma): with alpha
#' normalized to 1, the steady state S = ((1-theta)/theta) I + (1/gamma)
#' (2 - 1/theta) (Eq 14) is fitted by weighted least squares per TAD.
#' One-sided parameterization phi(alpha, beta, gamma): the steady state
#' S = (gamma/beta) I (Eq 13) identifies only the ratio gamma/beta,
#' reported as \code{$ratio}.
#'
#' @param S,I TAD x cell matrices of contacts in the ~50 kb windows of the
#'   S (stabilizing) and I (initiation) TAD boundaries
#' @param mode "two-sided", "one-sided", or "auto" (per-TAD modes taken
#'   from an \code{lv_si} object supplied as \code{mode_tbl})
#' @param delta_t projection horizon (arbitrary time units)
#' @param kCells cell pooling neighbourhood size (1 = no pooling)
#' @param cellKNN optional precomputed pooling graph
#' @param fit.quantile keep the top/bottom fraction of cells in the fit
#'   (NULL = all cells)
#' @param diagonal.quantiles weight on the S+I diagonal instead of the I
#'   marginal
#' @param mode_tbl output of [assign_si_boundaries] for mode = "auto"
#' @param mult size normalization target (contacts per \code{mult})
#' @param n.cores worker processes
#' @param verbose progress messages
#' @return list with \code{$theta}, \code{$gamma} (two-sided),
#'   \code{$ratio} (one-sided), \code{$mode}, \code{$ko}, normalized
#'   \code{$S}, \code{$I}, \code{$current}, \code{$projected},
#'   \code{$deltaE}, \code{$cellKNN}, \code{$mult}, \code{$delta_t}
#' @examples
#' # cells on the two-sided steady-state constraint line S = ((1-th)/th) I + b
#' set.seed(1); n <- 60
#' th <- 0.7; Ig <- runif(n, 0.5, 3)
#' S <- matrix((1 - th) / th * Ig + 2, 4, n,
#'            dimnames = list(paste0("tad", 1:4), paste0("c", 1:n)))
#' I <- matrix(Ig, 4, n, dimnames = dimnames(S))
#' # mult = 1e9 makes both size factors equal and constant, so the fit
#' # recovers the constraint slope exactly (alpha = 1 convention)
#' vel <- loop_velocity_estimates(S, I, mode = "two-sided", kCells = 1,
#'                                mult = 1e9)
#' median(vel$theta)
#' @export
loop_velocity_estimates <- function(S, I, mode = c("two-sided", "one-sided", "auto"),
                                    delta_t = 1, kCells = 10, cellKNN = NULL,
                                    fit.quantile = NULL, diagonal.quantiles = FALSE,
                                    mode_tbl = NULL, mult = 1e3,
                                    n.cores = default_n_cores(), verbose = TRUE) {
  mode <- match.arg(mode)
  if (!all(colnames(S) == colnames(I))) stop("S and I must share columns (cells)")
  vt <- intersect(rownames(S), rownames(I))
  if (length(vt) == 0) stop("S and I share no TADs (rows)")
  S <- S[vt, , drop = FALSE]; I <- I[vt, , drop = FALSE]

  # --- per-cell size normalization, separate factors per matrix ---
  # Each boundary matrix is normalized by its own cell totals (contacts
  # per `mult`), pooled over the cell kNN graph, then divided by the
  # pooled factors. This follows the velocity-framework preprocessing
  # the supplementary notes reference.
  S.cs <- pmax(Matrix::colSums(S) / mult, 1e-3)
  I.cs <- pmax(Matrix::colSums(I) / mult, 1e-3)
  if (kCells > 1) {
    if (is.null(cellKNN)) {
      if (verbose) cat("building cell kNN pooling graph ... ")
      ln <- t(log(t(t(I) / I.cs) + 1))   # cells x TAD-features
      # n.threads = 1: an OpenMP thread pool must not be active before the
      # mclapply forks below (forked children inherit locked mutexes and
      # silently corrupt results); per-TAD parallelism covers the cost
      cellKNN <- balanced_knn(ln, kCells, kCells * 100, n.threads = 1)
      diag(cellKNN) <- 1
      if (verbose) cat("done\n")
    }
    cm <- colnames(S)
    S.csp <- as.numeric(S.cs %*% cellKNN[cm, cm])
    I.csp <- as.numeric(I.cs %*% cellKNN[cm, cm])
    Sn <- t(t(as.matrix(S %*% cellKNN[cm, cm])) / S.csp)
    In <- t(t(as.matrix(I %*% cellKNN[cm, cm])) / I.csp)
  } else {
    cellKNN <- NULL
    Sn <- t(t(as.matrix(S)) / S.cs); In <- t(t(as.matrix(I)) / I.cs)
  }
  I.csn <- if (kCells > 1) I.csp else I.cs   # for the offset threshold

  if (mode == "auto") {
    if (is.null(mode_tbl)) {
      warning("mode = 'auto' without mode_tbl; defaulting to two-sided")
      modes <- rep("two-sided", length(vt))
    } else {
      modes <- mode_tbl$mode[match(vt, rownames(mode_tbl$S))]
      modes[is.na(modes)] <- "two-sided"
    }
  } else modes <- rep(mode, length(vt))
  names(modes) <- vt

  if (verbose) cat("fitting steady-state models per TAD ... ")
  fits <- parallel::mclapply(seq_along(vt), function(k) {
    s <- as.numeric(Sn[k, ]); i <- as.numeric(In[k, ])
    if (!is.null(fit.quantile)) {
      # quantile-weighted intercept fit; the intercept doubles as the
      # offset subtracted from S before projection (paper pipeline)
      w <- rep(1, length(s))
      if (diagonal.quantiles) {
        smax <- max(stats::quantile(s, 0.99), 1e-3)
        imax <- max(stats::quantile(i, 0.99), 1e-3)
        x <- s / smax + i / imax
        eq <- stats::quantile(x, c(fit.quantile, 1 - fit.quantile))
        w <- as.numeric(x >= eq[2] | x <= eq[1])
      } else {
        eq <- stats::quantile(i, c(fit.quantile, 1 - fit.quantile))
        w <- as.numeric(i >= eq[2] | i <= eq[1])
      }
      if (sum(w > 0) < 5) w <- rep(1, length(s))
      fl <- lv_wls_fit(i, s, w)   # S = a*I + b (Eq 14: S as the response)
      o <- fl$intercept; a <- fl$slope
    } else {
      # unweighted-by-quantile path: offset from cells with near-zero I
      # (threshold one raw count in normalized units), then a
      # fourth-moment-weighted through-origin slope on the offset-adjusted S
      o <- 0
      zi <- i < 1 / I.csn
      if (any(zi)) o <- sum(s[zi]) / (sum(zi) + 1)
      w <- i^4 + s^4
      y <- s - o
      den <- sum(w * i * i)
      a <- if (is.finite(den) && den > 0) sum(w * i * y) / den else NA_real_
    }
    if (!is.finite(a))
      return(c(theta = NA_real_, gamma = NA_real_, ratio = NA_real_,
               geff = NA_real_, o = o))
    if (modes[k] == "two-sided") {
      theta <- 1 / (1 + max(a, 1e-3))
      theta <- min(max(theta, 0.51), 0.99)   # identifiability clamp
      gamma <- if (is.finite(o) && o > 1e-8) (2 - 1 / theta) / o else NA_real_
      c(theta = theta, gamma = gamma, ratio = NA_real_, geff = a, o = o)
    } else {
      # Eq 13: S = (gamma/beta) I -> the slope itself is the ratio
      c(theta = NA_real_, gamma = NA_real_, ratio = max(a, 1e-3),
        geff = a, o = o)
    }
  }, mc.cores = n.cores, mc.preschedule = TRUE)
  ko <- do.call(rbind, fits); rownames(ko) <- vt
  if (verbose) cat("done\n")

  # Steady-state extrapolation (paper form): the per-cell velocity signal is
  # the offset-corrected residual of the observed state from the fitted
  # steady-state constraint, saturated over the horizon,
  #   dI = (1 - exp(-geff * delta_t)) * (max(S - o, 0) / geff - I)
  # with geff the unclamped S ~ I regression slope (alpha is not
  # identifiable in the steady state, so the observed S level acts as the
  # generation scale) and o the per-TAD basal-contact offset.
  valid <- is.finite(ko[, "geff"])
  vv <- vt[valid]
  deltaE <- t(vapply(vv, function(k) {
    geff <- ko[k, "geff"]
    y <- pmax(Sn[k, ] - ko[k, "o"], 0)
    as.numeric((1 - exp(-geff * delta_t)) * (y / geff - In[k, ]))
  }, numeric(ncol(In))))
  rownames(deltaE) <- vv; colnames(deltaE) <- colnames(In)

  # the observed I values represent the current cell state and the
  # extrapolated I values the projected state (supplementary notes,
  # "Inference of parameters and extrapolation of cell state"); the current
  # state uses the unpooled size-normalized matrix restricted to the
  # fitted TADs, as in the velocity framework the notes reference
  Icur <- t(t(as.matrix(I[vv, , drop = FALSE])) / I.cs)
  proj <- Icur + deltaE; proj[proj < 0] <- 0
  list(theta = ko[, "theta"], gamma = ko[, "gamma"], ratio = ko[, "ratio"],
       mode = modes, ko = ko, S = Sn, I = In,
       current = Icur, projected = proj, deltaE = deltaE,
       cellKNN = cellKNN, mult = mult, delta_t = delta_t)
}
