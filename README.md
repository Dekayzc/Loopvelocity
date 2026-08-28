# Loopvelocity

Loop velocity estimation from single-cell Hi-C.

Loop velocity models TAD-boundary contact dynamics with a loop-extrusion
kinetic model — one-sided `phi(alpha, beta, gamma)` or two-sided
`phi(alpha, theta, gamma)` — and extrapolates cell chromatin states, in the
spirit of RNA velocity but for contacts.

## Install

    remotes::install_github("Dekayzc/Loopvelocity")

## Quick start

    library(Loopvelocity)
    si <- read_tad_contacts("tad_left_50k.bed", "tad_right_50k.bed")
    vel <- loop_velocity_estimates(si$S, si$I, mode = si$mode, kCells = 10,
                                   fit.quantile = 0.05)
    vel <- loop_velocity_em(vel)
    vel <- project_cells(vel, delta_t = 1)
    plot_loop_velocity(umap_embedding, vel, kernel = "correlation",
                       show.grid.flow = TRUE)

## Full Hi-C pipeline

    cells <- read_scool("cells.scool")
    bnd   <- call_tad_boundaries(contact_matrix, window = 5)
    si    <- assign_si_boundaries(boundary_contacts(cells, bnd))

## Models

- **Steady state**: per-TAD weighted linear fit of the constraint lines
  (Eq 13/14 of the supplementary notes); two-sided estimates `theta` and
  `gamma` under `alpha = 1`, one-sided the ratio `gamma/beta`.
- **Dynamical (EM)**: per-cell time as hidden variable, E-step and M-step
  both solved with box-constrained L-BFGS (Eq 17–19); initialized from the
  steady-state fit.
- **Projection**: `dI = (1/gamma - I) (1 - exp(-(1-theta) gamma dt))`
  (Eq 16, two-sided) or forward evaluation of the one-sided analytic
  solution (Eq 5–6).
