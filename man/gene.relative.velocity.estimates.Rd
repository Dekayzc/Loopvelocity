% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/momentum_routines.R
\name{gene.relative.velocity.estimates}
\alias{gene.relative.velocity.estimates}
\title{Estimate RNA velocity using gene-relative slopes}
\usage{
gene.relative.velocity.estimates(
  emat,
  imat,
  deltaT = 1,
  smat = NULL,
  steady.state.cells = colnames(emat),
  kCells = 10,
  cellKNN = NULL,
  kGenes = 1,
  theta.s = FALSE,
  RNA_mat = NULL,
  old.fit = NULL,
  mult = 1000,
  min.imat.smat.correlation = 0.05,
  min.imat.emat.correlation = 0.05,
  min.imat.emat.slope = 0.05,
  zero.offset = FALSE,
  deltaT2 = 1,
  fit.quantile = NULL,
  diagonal.quantiles = FALSE,
  show.gene = NULL,
  do.par = TRUE,
  cell.dist = NULL,
  emat.size = NULL,
  imat.size = NULL,
  cell.emb = NULL,
  cell.colors = NULL,
  expression.gradient = NULL,
  residual.gradient = NULL,
  n.cores = defaultNCores(),
  verbose = TRUE
)
}
\arguments{
\item{emat}{- exons count matrix}

\item{imat}{- introns count matrix}

\item{deltaT}{- amount of time to project the cell forward}

\item{smat}{- optional spanning read matrix (used in offset calculations)}

\item{steady.state.cells}{- optional set of steady-state cells on which the gamma should be estimated (defaults to all cells)}

\item{kCells}{- number of k nearest neighbors (NN) to use in slope calculation smoothing}

\item{cellKNN}{- optional pre-calculated cell KNN matrix}

\item{kGenes}{- number of genes (k) to use in gene kNN pooling}

\item{theta.s}{- option to calculate value of theta}

\item{RNA_mat}{- RNA count matrix for theta calculation}

\item{old.fit}{- optional old result (in this case the slopes and offsets won't be recalculated, and the same kNN graphs will be used)}

\item{mult}{- library scaling factor (1e6 in case of FPM)}

\item{min.imat.smat.correlation}{- minimum required Spearman rank correlation between n and s counts of a gene}

\item{min.imat.emat.correlation}{- minimum required Spearman rank correlation between n and e counts of a gene}

\item{min.imat.emat.slope}{- minimum sloope of n~e regression}

\item{zero.offset}{- should offset be set to zero, or determined (through smat regression or using near-0 e cases)}

\item{deltaT2}{- scaling of the projected difference vector (normally should be set to 1)}

\item{fit.quantile}{perform gamma fit on a top/bottom quantiles of expression magnitudes}

\item{diagonal.quantiles}{whether extreme quantiles should be computed diagonally}

\item{show.gene}{an optional name of a gene for which the velocity estimation details should be shown (instead of estimating all velocities)}

\item{do.par}{whether the graphical device parameters should be reset as part of show.gene (default=TRUE)}

\item{cell.dist}{- cell distance to use in cell kNN pooling calculations}

\item{emat.size}{- pre-calculated cell sizes for the emat (spliced) matrix}

\item{imat.size}{- pre-calculated cell sizes for the imat (unspliced) matrix}

\item{cell.emb}{- cell embedding to be used in show.gene function}

\item{cell.colors}{- cell colors to be used in show.gene function}

\item{expression.gradient}{- color palette used to show the expression magnitudes in show.gene function}

\item{residual.gradient}{- color palette used to show the u residuals in show.gene function}

\item{n.cores}{- number of cores to use}

\item{verbose}{- output messages about progress}
}
\value{
a list with velocity results, including the current normalized expression state ($current), projected ($projected) over a certain time ($deltaT), unscaled transcriptional change ($deltaE), fit results ($gamma, $ko, $sfit if spanning reads were used), optional cell pooling parameters ($cellKNN, $kCells), kNN-convolved normalized matrices (conv.imat.norm and conv.emat.norm), library scale ($mult)
}
\description{
Estimate RNA velocity using gene-relative slopes
}
\examples{
\dontrun{
 # use min/max quantile gamma fit (recommended option when one can afford to do cell kNN smoothing)
 # The example below uses k=5 cell kNN pooling, and top/bottom 2\% exprssion quantiles
 # emat and imat are spliced (exonic) and unspliced (intronic) molecule/read count matirces
(preferably filtered for informative genes)
 rvel <- gene.relative.velocity.estimates(emat,imat,deltaT=1,kCells = 5,fit.quantile = 0.02)

 # alternativly, the function can be used to visualize gamma fit and regression for a
particular gene. here we pass embedding (a matrix/data frame with rows named with cell names,
and columns corresponding to the x/y coordinates)

 # and cell colors. old.fit is used to save calculation time.
 gene.relative.velocity.estimates(emat,imat,deltaT=1,kCells = 5,fit.quantile = 0.02,
    old.fit=rvel,show.gene='Chga',cell.emb=emb,cell.colors=cell.colors)
}
}
