# Loopvelocity

**Loop velocity: estimating the time derivative of the chromatin contact
state from single-cell Hi-C.**

Loop velocity adapts the idea of RNA velocity to the spatial organization
of the genome. Instead of following transcriptional kinetics, it models the
kinetics of **loop extrusion**: contacts around TAD boundaries drift
through time as extrusion complexes load, translocate and dissociate, and
this drift indicates where a cell's chromatin contact state is heading —
its future on a developmental trajectory.

Given per-cell contact counts in ~50 kb windows around the two boundaries
of each TAD, Loopvelocity fits a loop-extrusion kinetic model, extrapolates
each cell's contact state forward in time, and visualizes the resulting
flow on any cell embedding (UMAP / tSNE / PCA).

## Method

Contacts at the **S (stabilizing)** and **I (initiation)** TAD boundaries
follow two alternative kinetic schemes, chosen per TAD by a significance
test of the two boundaries' contact numbers (a significant difference
indicates two-sided extrusion; the higher-contact boundary is S).

**One-sided extrusion** — the S boundary stays put and actively transfers
contacts toward the I boundary at rate β, contacts form at rate α, and the
loop dissociates at rate γ:

    dS/dt = α − βS        S(t) = (α/β)(1 − e^(−βt))
    dI/dt = βS − γI       I(t) = (α/γ)(1 − e^(−γt)) + (α/(γ−β))(e^(−γt) − e^(−βt))

Parameters: **φ(α, β, γ)**.

**Two-sided extrusion** — DNA translocates through the extrusion complex
at both boundaries, so there is no directed S→I contact transfer; instead
the shared formation rate α is split between the boundaries by **θ** and
both boundaries share the dissociation rate γ:

    dS/dt = θα − θγS      S(t) = (α/γ)(1 − e^(−θγt))
    dI/dt = (1−θ)α − (1−θ)γI    I(t) = (α/γ)(1 − e^(−(1−θ)γt))

Parameters: **φ(α, θ, γ)**.

**Estimation.** The steady state model fits the equilibrium constraint of
each scheme by weighted linear regression per TAD (balanced-kNN cell
pooling and top/bottom-quantile weighting follow the velocity framework
referenced in the paper): two-sided gives θ from the slope and γ from the
intercept under the α = 1 normalization; one-sided identifies the ratio
γ/β. The dynamical model then refines the parameters with an EM algorithm
over the full analytic solutions: per-cell time is the hidden variable,
the E-step (hidden times) and M-step (parameters) are both solved with
box-constrained L-BFGS, and the rate scale is fixed by the α = 1
convention (absolute rates are otherwise unidentifiable because
rescaling all rates and time jointly leaves the trajectories unchanged).

**Extrapolation.** The projected cell state is the extrapolated I:
ΔI = (1 − e^(−γ_eff·Δt)) · (max(S − o, 0)/γ_eff − I), where γ_eff is the
per-TAD S ~ I regression slope and o the per-TAD basal-contact offset —
the velocity signal of each cell is its residual from the steady-state
constraint.

## Install

```r
# release
remotes::install_github("Dekayzc/Loopvelocity")

# development
remotes::install_github("Dekayzc/Loopvelocity", ref = "main")
```

System requirements: a C++14 compiler; OpenMP for multithreaded kernels.
Optional: `hdf5r` (cooler-format input), `strawr` (Juicer .hic input),
`Rtsne` (tSNE diagnostic plot).

## Quick start (pre-summarized boundary tables)

If you already have per-cell contact counts around TAD boundaries as two
tables (chr, start, end, contacts, cell, regionid — rows pair across the
two files by chromosome and region id):

```r
library(Loopvelocity)

# S/I assignment: higher-contact side becomes S; Wilcoxon test picks the
# per-TAD extrusion mode (use test = "t" for the t-test variant)
si <- read_tad_contacts("tad_left_50k.bed", "tad_right_50k.bed",
                        test = "wilcox", padj = 0.05)
sig <- si$mode == "two-sided"          # keep significant TADs

# steady-state loop velocity (what the paper's Figure 4 flow fields use)
vel <- loop_velocity_estimates(si$S[sig, ], si$I[sig, ], mode = "two-sided",
                               kCells = 10, fit.quantile = 0.05, n.cores = 8)

# flow field on a cells x 2 embedding (rownames = cell ids)
plot_loop_velocity(umap_embedding, vel, kernel = "correlation",
                   show.grid.flow = TRUE, grid.n = 40, arrow.scale = 10)
```

Optional dynamical refinement:

```r
vel <- loop_velocity_em(vel, iteration = 100, tol = 0.01, n.cores = 8)
vel <- project_cells(vel, delta_t = 1)
```

## Full Hi-C pipeline (raw contacts to flow field)

```r
library(Loopvelocity)

# 1. read cells (4DN single-cell cooler; also read_cool/read_mcool/
#    read_pairs/read_sparse_contacts/read_hic)
cells <- read_scool("cells.scool")

# 2. call TAD boundaries from an insulation profile
bnd <- call_tad_boundaries(contact_matrix, window = 5, min_delta = 0.1)

# 3. per-cell contacts in the ~50 kb windows around each boundary
bc <- boundary_contacts(cells, bnd, flank = 25e3)

# 4. S/I roles + extrusion mode per TAD
si <- assign_si_boundaries(bc$left, bc$right, test = "t", padj = 0.05)

# 5. steady state -> EM -> projection -> visualization
vel <- loop_velocity_estimates(si$S, si$I, mode = si$mode,
                               kCells = 10, fit.quantile = 0.05)
vel <- loop_velocity_em(vel)
vel <- project_cells(vel, delta_t = 1)
plot_loop_velocity(emb, vel, show.grid.flow = TRUE)
```

## Function reference

**Input & preprocessing**

| Function | Purpose |
|---|---|
| `read_tad_contacts()` | boundary BED tables → S/I matrices + mode test (quick path) |
| `read_scool()` / `read_mcool()` / `read_cool()` | cooler-family contact files (via hdf5r) |
| `read_pairs()` | 4DN .pairs contact lists |
| `read_sparse_contacts()` | HiC-Pro/Homer sparse triplet matrices |
| `read_hic()` | Juicer .hic pixels (via strawr) |
| `import_boundaries()` | external TAD boundary tables |
| `compute_insulation_score()` / `call_tad_boundaries()` | insulation-based boundary calling |
| `boundary_contacts()` | per-cell windowed contact counts around boundaries |
| `assign_si_boundaries()` | S/I role assignment + one/two-sided mode test |
| `filter_tads_by_activity()` | drop inactive TADs |

**Estimation**

| Function | Purpose |
|---|---|
| `loop_velocity_estimates()` | steady-state fits (Eq 13/14) with kNN pooling and quantile weighting |
| `loop_velocity_em()` | EM dynamical model (hidden cell time; L-BFGS E/M steps) |
| `project_cells()` | extrapolate cell states (steady-state or EM parameters) |
| `lv_two_sided_solution()` / `lv_one_sided_solution()` | evaluate the analytic ODE solutions |
| `lv_project_two_sided()` | two-sided projection increment (Eq 16) |

**Visualization**

| Function | Purpose |
|---|---|
| `plot_loop_velocity()` | flow field on an embedding; correlation or euclidean kernel |
| `plot_velocity_pcs()` | diagnostic: consecutive PC pairs with velocity arrows |
| `plot_velocity_tsne()` | tSNE view of the velocity field (Rtsne) |
| `plot_si_fit()` | per-TAD S–I scatter with the fitted constraint line |

## Output structure

`loop_velocity_estimates()` returns a list:

- `$theta`, `$gamma` — two-sided steady-state parameters (α = 1 scale);
  `$ratio` — one-sided γ/β; NA where a TAD's fit is degenerate
- `$mode` — per-TAD `"two-sided"` / `"one-sided"`
- `$ko` — per-TAD coefficients (theta, gamma, ratio, geff, offset o)
- `$S`, `$I` — kNN-pooled size-normalized boundary matrices
- `$current` — the observed I state (unpooled, size-normalized)
- `$deltaE` — the velocity signal (projected increment per TAD per cell)
- `$projected` — extrapolated cell state
- `$cellKNN` — the pooling graph (reuse via `cellKNN =` to save time)

`loop_velocity_em()` adds `$par` (per-TAD θ, γ or β, γ), `$time` (hidden
per-cell times), `$iterations`, `$converged`.

## Testing

The package ships a testthat suite covering the analytic solutions
(against an independent numerical integration), parameter recovery from
simulated trajectories, all readers (on synthetic files), and plotting
smoke tests: run `devtools::test()` or `R CMD check`.

## Citation

Zhang, Chen, et al. "Exploring DNA movement through the application of
droplet based high efficient chromatin conformation capture (DropHiChew)
and loop velocity." *bioRxiv* (2024): 2024-06.

## License

GPL-3.
