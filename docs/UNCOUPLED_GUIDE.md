# Atlas Pipeline — `run_atlas_surface_pipeline.m`

Master reference for the HSS/HPSS/MPSS atlas pipeline: what it does, how to run
it, and how each function in it works. For the models and their parameters see
[MODELS](MODELS.md); the load-bearing modelling rules are in section 10. In this
package the script sits in the package folder; its outputs go to `data/`
(git-ignored). The coupled pipeline's wrapper builds only Stages 1-3 (the seed
libraries) through `atlas_build_libraries`.

---

## 1. What the pipeline produces

From one reaction–diffusion parameter set it builds a **nine-sheet atlas** of
half-domain steady states on the mass-partition grid `(Δn_PAK, Δn_GEF)`, then
runs one PDE trajectory and measures how well the atlas explains it.

The model has five abstract components `(u, q, v, w, p)`: two slow-diffusing
components and three fast-diffusing components. Each half-domain steady state
is one
of three **branch types**:

| Type | Meaning | How it is obtained |
|------|---------|--------------------|
| **HSS**  | Homogeneous steady state (spatially uniform) | algebraic multi-root solve |
| **HPSS** | Half-peak steady state (peak at the domain boundary) | 14-variable BVP on `[0, L/2]` |
| **MPSS** | Mid-peak steady state (interior peak) | **constructed**: two quarter-domain HPSS solves reflected peak-to-peak — never a direct interior-peak BVP |

A **sheet** `<L>+<R>` stitches a left half of type `L` to a right half of type
`R` (mirror cell), giving 9 sheets: HSS_HSS, HSS_HPSS, HSS_MPSS, HPSS_HSS,
HPSS_HPSS, HPSS_MPSS, MPSS_HSS, MPSS_HPSS, MPSS_MPSS. Each sheet carries its own
redistribution-flux field and nullclines on the shared grid.

---

## 2. Standard usage

```matlab
% edit the parameter + settings block at the top of the file, then:
run_atlas_surface_pipeline
```

Everything is driven by the settings block at the top of the script (§4). A run
executes Stages 1–12 in order; each stage either rebuilds its output or reuses
the existing file (§5). Every run writes into its own per-parameter-set folder
`data/<params>/` (`outputs/` + `plots/` subfolders); the trajectory-specific
Stage 8–12 outputs additionally go into a per-initial-condition subfolder
`data/<params>/<sheet>_<i>_<j>/` so different ICs at the same parameters never
overwrite each other (§8). Figures are written regardless of whether windows
are shown on screen.

**Typical workflows**

- *First run / new parameter set:* set the four `force_*` flags to `true` so the
  branch libraries and sheets rebuild from scratch.
- *Re-examine an existing atlas with a different PDE trajectory:* leave `force_*`
  `false`, edit the `trace_*` settings — Stages 8–9 auto-detect the change and
  rerun only the trajectory + classification (§5).
- *Re-style plots only:* edit `trace_n_top_sheets` / `movie_*` / `show_figures`
  — these are display-only and never trigger a resolve.

---

## 3. Stage sequence

| # | Stage | Function | Output |
|---|-------|----------|--------|
| 1 | HPSS branch | `run_hybrid_sweeper_pdefb` → `build_branch_library` | `bvp_sweep_full_results.mat`, `branch_HPSS_latest.mat` |
| 2 | HSS branch | `sweep_hss_multibranch_new` → `build_branch_library` | `hss_sweep_multibranch_results.mat`, `branch_HSS_latest.mat` |
| 2b | HSS dispersion stability | `plot_hss_dispersion_atlas` | `HSS_Dispersion_Stability.pdf`, `hss_dispersion_atlas.mat` |
| 3 | MPSS branch | `build_mpss_branch` | `branch_MPSS_latest.mat` |
| 4 | Branch surfaces | `visualize_atlas_surfaces` | branch surface PDFs |
| 5 | Stitch sheets | `build_atlas_sheets` | `atlas_sheets_latest.mat` |
| 6 | Per-sheet flux plots | `plot_atlas_sheets` | `Atlas_Sheet_*.pdf`, overview, combined nullclines |
| 7 | Combined 3-D surfaces | `plot_atlas_surfaces_combined` | `Atlas_Surfaces_Combined.pdf` |
| 8 | PDE trajectory | `run_atlas_pde_trace` | `atlas_pde_trace_latest.mat` |
| 9 | Classify trajectory | `classify_pde_to_atlas_sheets` | `pde_atlas_classification_latest.mat` |
| 10 | Phase-plane pages | `plot_atlas_phase_plane` | `Atlas_PhasePlane_*.pdf` |
| 11 | Classification summary | `plot_pde_to_atlas_classification` | `PDE_Atlas_Classification.pdf`, `PDE_Atlas_Frame_Overlays.pdf` |
| 12 | Trajectory movie | `plot_atlas_pde_movie` | `Atlas_PDE_Movie.mp4` |

Stages 1–3 build the branch libraries; 4–7 visualize the static atlas; 8–12
run and analyse one PDE trajectory against it. Stage 2b (display-only, runs
right after the HSS branch is built) is the linear-stability companion to the
HSS branch — see §6.

---

## 4. Settings reference (top of the script)

**Model parameters** — `p.a … p.f`, `p.Dm`, `p.Dc`, `p.L`, `p.M_BG`, `p.M_P`,
`p.res`. The reaction-rate and diffusion constants, domain length, total GEF and
PAK masses, and grid resolution. The shipped script holds the reference parameter
set. `p.res` is the grid side length (e.g. 50 → 2500 cells).

**Solver options** (only consulted when a sweep/branch actually rebuilds):
- `p.seeding_mode` — `'auto'` (headless) or `'interactive'` (pick the seed
  visually).
- `p.skip_stability` — `true` skips the per-cell eigenvalue stability check in
  the HPSS sweep (leaves stability uncomputed). See §7.
- `p.mpss_strategy` — MPSS construction strategy (`'continuation'`).
- `p.phase0_construct` (default `true`) — HPSS sweep (Stage 1) Phase 0 Tier 0:
  seed by solving the BVP directly from a constructed analytic edge-peak
  guess, no PDE relaxation. Reaches a dynamically-unstable static half-peak
  that relaxation-based seeding cannot (relaxation always flees an unstable
  steady state). Falls through to the ordinary relaxation/L-creep seed if no
  constructed guess validates.
- `p.osc_fallback` (default `false`) — oscillation-aware PDE-relaxation
  fallback for Stage 1's Phase 2 and Stage 3's quarter-domain continuation:
  detects a non-relaxing (limit-cycle) chunk sequence and offers time-mean /
  phase-snapshot `bvp4c` seed candidates instead of an arbitrary limit-cycle
  frame. Needed for parameter sets that drive sustained PDE oscillation. Key
  tuning knob is `p.pde_T_chunk` (chunk length must span several oscillation
  periods to make non-relaxation detectable) — measure the period with a
  short standalone `pdepe` probe before setting it for a new parameter set.
  The development repository holds the mechanism notes and validated numbers.

**Rebuild flags** — `force_hpss`, `force_hss`, `force_mpss`, `force_sheets`.
Each forces its stage (1, 2, 3, 5) to rebuild even if the output file exists.
Leave `false` to reuse existing files. `force_trace` forces Stages 8–9; normally
unneeded because they auto-detect staleness (§5).

**Trajectory settings** (Stage 8, part of the staleness key):
- `trace_sheet`, `trace_cell` — which sheet and grid cell `[i_GEF, j_PAK]` to
  seed the PDE initial condition from.
- `trace_T_end`, `trace_t_fine_end`, `trace_dt_fine`, `trace_dt_coarse` — a
  **two-phase linear time grid**: step `trace_dt_fine` over `[0, trace_t_fine_end]`
  to resolve the fast early transient, then the coarser `trace_dt_coarse` out to
  `trace_T_end`. Assembled into `trace_t_span`.
- `trace_N_mesh` — spatial mesh points over `[0, L]`.

**Classification thresholds** (Stage 9, part of the staleness key):
- `trace_ambiguity_gap`, `trace_off_atlas_tol` — provisional residual thresholds
  for flagging a frame ambiguous or off-atlas. Not calibrated decisions; the
  continuous residuals are the primary output.

**Display-only settings** (never trigger a resolve):
- `trace_n_top_sheets` (Stage 10) — `[]` draws one page per dwelled sheet; `N`
  draws the N best-fitting sheets by min residual (even non-matching ones);
  `Inf`/9 draws all nine.
- `render_movie`, `movie_fps`, `movie_frame_stride`, `movie_resolution_scale`,
  `movie_quality` (Stage 12).
- `show_figures` — `true` pops figure windows on screen; `false` is fully
  headless. Either way every figure is exported to `data/<params>/plots/`.

---

## 5. Rebuild / caching model

**Stages 1–3, 5** rebuild only when their `force_*` flag is `true` or their
output `.mat` is missing; otherwise they reuse the existing file. Rebuilding a
branch changes its `computation_id`.

**Stages 8–9** are settings-aware. Stage 8 packs the trajectory-defining settings
(`trace_sheet`, `trace_cell`, `trace_t_span`, `trace_N_mesh`, and the source
branch/atlas `computation_id`s) into `trace_cfg` and saves it inside the trace
`.mat` as `run_config`. On the next run `decide_trace` compares the current
`trace_cfg` to the stored one; any difference (including a rebuilt upstream
branch) regenerates the trace. Stage 9 (`decide_class`) reclassifies iff the
trace was regenerated, the classification is missing, it was built from a
different trace, or the thresholds changed.

**Stages 4, 6, 7, 10, 11, 12** always run — they only display already-computed
data and have no staleness key.

Because a branch rebuild changes its id, setting `force_hpss/hss/mpss` also
cascades a trace regeneration → reclassification → all downstream plots.

---

## 6. Function reference

### Stage 1 — `run_hybrid_sweeper_pdefb(p)`
Solves the 14-variable HPSS BVP at every grid cell, sweeping mass partition at
fixed `L`. The centre seed is found first by a **constructed-guess cold start**
(`p.phase0_construct`, default on — solves the BVP directly from an analytic
edge-peak guess, no relaxation), falling through to relaxation/L-creep seeding
if no guess validates; the seed is then continued across the grid by
predictor–corrector, and where that fails it uses a **PDE-relaxation fallback**
(extrapolate a neighbour profile, relax it under the reaction–diffusion PDE —
oscillation-aware if `p.osc_fallback` — feed the relaxed state to `bvp4c`).
Failed cells are masked, never fabricated. Phase 5 (optional, gated by
`p.skip_stability`) classifies each solved cell's temporal stability via a
discretised-Jacobian eigenvalue check.

Every converged candidate passes the shared morphology classifier
`classify_hpss_morphology` before it can be stored or reused as a continuation
seed. A valid HPSS may have a slight boundary dip and a crown within the first
10% of the half-domain (a half-volcano), but the opposite boundary must remain
the trough and there may be no additional significant crown or interior trough.
Before Phase 4 and after every recovery pass, the complete seed pool is audited;
HSS collapses, interior-peak/MPSS-like profiles, multi-peak profiles, and
unclassified profiles are quarantined as failed HPSS cells. Thus recovery may
approach the HPSS existence boundary but cannot propagate a different branch
through it.

Output `bvp_sweep_full_results.mat`: `PAK_DIFF`, `GEF_DIFF`, `{U,Q,V,W,P}_trough`,
`SOL_grid` (cell of BVP structs), `failure_mask`, `is_stable_grid`, geometry,
`par`, `p`, `computation_id`, plus `morphology_class_grid`,
`morphology_reason_grid`, `morphology_diagnostics_grid`, and
`quarantined_nonhpss_mask`. The quarantine mask records cells where the
grid-level or Phase-4 safety gate encountered a converged non-HPSS candidate; it
does not assert that an HPSS branch is mathematically nonexistent there.

### Stage 2 — `sweep_hss_multibranch_new(p)`
Finds all real-positive homogeneous steady states per cell by multistart
`fsolve`, then assigns **global branch labels** by nearest-neighbour continuity
in the full `[u q v w p]` state (not by per-cell root order), so a physical
branch keeps one label across the grid. Output
`hss_sweep_multibranch_results.mat`: `{U,Q,V,W,P}_hss_branches` (`res×res×maxK`),
`n_branches_grid`, `branch_id_grid`, `branch_summary`, provenance.

### Stages 1–2 adapter — `build_branch_library(hpss_path, hss_path, out_dir)`
Read-only adapter mapping the HPSS and HSS sweep files into the common
**branch-library schema** used by everything downstream. HPSS existence is
morphology-first: a cell counts as HPSS only if the solver converged **and** the
profile is a genuine boundary half-peak (amplitude ≥ 0.5), so HSS-collapsed
profiles cannot contaminate an HPSS sheet. It uses the same shared classifier as
the sweeper and emits per-cell `branch_class` and `morphology_reason` alongside
`existence_mask`, `SOL`, `flux_values`, `flux_source_location`, `amp`, `peak_x`,
`mass_residual`, and annotation-only stability fields. Writes
`branch_HPSS_latest.mat` and `branch_HSS_latest.mat`; never modifies the source
files. (Called after Stage 1 and again after Stage 2.)

### Stage 2b — `plot_hss_dispersion_atlas(opts)`
Linear-stability atlas of the per-cell HSS by **dispersion relation**. At each
cell it forms the reaction Jacobian `J` of the representative HSS branch
(analytic; verified to satisfy both mass-conservation null-vectors
`[1 1 1 1 0]`, `[0 1 0 0 1]` exactly), scans the discrete no-flux modes
`kₙ = nπ/ℓ` (`ℓ = L/2`, half-domain; `n ≥ 1` only, which respects mass
conservation and sidesteps the two exactly-marginal `k=0` modes), and records
where the HSS is Turing-unstable (`max Re eig(J − kₙ²D) > tol`) and the
fastest-growing **onset** mode `n*`. Two panels over `(Δn_PAK, Δn_GEF)`:
(A) the unstable region shaded with opacity ∝ onset, the colorbar in **spatial
periods on the half-domain = `n*/2`** (0.5 = one half-peak = HPSS morphology,
1.0 = one full period, …; wavelength `Λ = L/n*`, so the reflected full-domain
count is `n*` peaks); (B) onset **character** — stationary Turing vs oscillatory
(Hopf–Turing → the shuttle). **`n*/2` is the linear ONSET, not the coarsened
attractor peak count** (mass-conserved systems coarsen — Invariant 9 spirit).
opts: `hss_path`, `save_dir` (required), `domain` (`'half'` default | `'full'`),
`n_max` (80), `tol` (1e-6), `visible`. Display-only; reads `branch_HSS_latest.mat`,
writes `HSS_Dispersion_Stability.pdf` + `hss_dispersion_atlas.mat` (grids
`unstable`, `nstar`, `onset_periods`, `maxgrowth`, `osc`, provenance).

### Stage 3 — `build_mpss_branch(p)`
Constructs the MPSS branch: a quarter-domain HPSS sweep, then peak-to-peak
reflection into a half-domain interior-peak profile. Flux is taken at the outer
trough on the half-domain interface. The quarter sweep's row continuation and
its multi-pass rescue both run as `parfeval`-per-task with a wall-clock timeout
(`p.phase3_row_timeout`, `p.cell_timeout`) plus checkpointing and resume
(Invariant 8) — a timed-out row's or cell's unsolved cells are simply left
failed for the rescue pass to pick up, nothing is fabricated. Failed quarter
solves or reflection checks are masked. Output `branch_MPSS_latest.mat` in the
same branch schema (MPSS stability is never computed — annotation fields are
all-false).

### Stage 4 — `visualize_atlas_surfaces(opts)`
Surface visualization of the three branch-library files: existence masks
(always raw), 3-D redistribution-flux surfaces (`v+w` and `p`, mirrored across the
origin), and scalar fields (`amp`, `peak_x`, `mass_residual`). Nullclines are
overlaid only if sheet files already exist. Display-only; reads branch files,
writes PDFs.

### Stage 5 — `build_atlas_sheets(hpss_path, hss_path, mpss_path, out_dir)`
Stitches the nine sheets by pure aggregation over the three branch files — no
new BVPs, no source file modified, no cell fabricated. The mirror of grid cell
`(i,j)` is exactly `(res+1-i, res+1-j)` (pure index arithmetic on the symmetric
grid, no interpolation). Per sheet: `sheet_exists = left.existence(i,j) &
right.existence(mirror)`; flux `U_flux = -2·Dc·(p_left − p_right)`,
`V_flux = -2·Dc·(vw_left − vw_right)` from each branch's own stored `flux_values`;
`sheet_stable_provisional = sheet_exists & half_stable_left & half_stable_right`
(annotation only). Output `atlas_sheets_latest.mat`: a `sheets` struct of nine
named sub-sheets plus `sheet_names`, geometry, `par`, and `source_*_id`
provenance.

### Stage 6 — `plot_atlas_sheets(opts)`
One phase-plane page per sheet — existence-region shading, normalized flux
quiver, and the two conserved-pool flux nullclines — reusing Stage 5's stored flux
(never recomputed). Also writes a 3×3 overview and a combined all-nullclines
figure.

### Stage 7 — `plot_atlas_surfaces_combined(opts)`
One 3-D figure overlaying all three branch flux surfaces with every sheet's
nullcline draped on top at its exact height (on a sheet's nullcline the left and
mirrored-right surface heights are equal by definition, so the height is
unambiguous). Reuses Stage 5's nullcline data.

### Stage 8 — `run_atlas_pde_trace(opts)`
Builds a full-domain PDE initial condition by reconstructing the chosen sheet's
stitched profile (via the shared `reconstruct_stitched_sheet`, so the t=0 frame
equals the classifier's reconstruction to machine precision), runs `pdepe` over
the given `t_span`, and records the per-frame mass-asymmetry trace
`(Δn_PAK(t), Δn_GEF(t))`. Output `atlas_pde_trace_latest.mat`: `sol_pde`
(`nt×nx×5`), `x_full`, `t_span`, `nPAK_trace`, `nGEF_trace`, `u0_stitched`,
IC/provenance fields, and the `run_config` staleness key.

### Stage 9 — `classify_pde_to_atlas_sheets(opts)`
For every PDE frame, snaps it to the nearest grid cell and computes, for each
sheet that exists there, the **relative-L2 spatial-profile residual** between the
frame and that sheet's reconstructed profile (Invariant 7 — profile residual,
not phase-plane distance). The best sheet is the argmin residual; frames are
flagged ambiguous (2nd-best within `ambiguity_gap`) or off-atlas (best residual
> `off_atlas_tol`). Output `pde_atlas_classification_latest.mat`: `t_span`,
`dpak`, `dgef`, `best_cell`, `residual_by_sheet` (`nt×9`),
`residual_by_species`, `best_sheet_idx/name`, `best/second_best_residual`,
`ambiguous`, `off_atlas`, and `provenance`.

### Stage 10 — `plot_atlas_phase_plane(opts)`
One page per selected sheet (dwelled sheets, or the top-`n_top` by min residual).
Each page shows that sheet's flux field/nullclines with the **full** trajectory
as a grey dotted background line and the on-this-sheet segments overlaid as a
thick 50%-opaque black line; the trace start is a grey dot, the end unmarked.
Requires the Stage 9 classification. Output `Atlas_PhasePlane_<sheet>.pdf`.

### Stage 11 — `plot_pde_to_atlas_classification(opts)`
The classifier's own diagnostic view: (1) best atlas sheet per frame over time
(confident/ambiguous/off-atlas colour-coded) and (2) the spatial-profile
residual of all nine sheets over time. Also writes a frame-vs-best-atlas-profile
overlay figure (all 5 species at a few times). Outputs
`PDE_Atlas_Classification.pdf` and `PDE_Atlas_Frame_Overlays.pdf`.

### Stage 12 — `plot_atlas_pde_movie(opts)`
Renders the raw PDE trajectory to `.mp4`: two stacked panels (`u,q` and
`v,w,p` vs `x`), one video frame per stored PDE frame (subsampled by
`frame_stride`), fixed per-species colours and fixed axis limits.
`resolution_scale` supersamples the raster (true vector video is not possible —
no standard container stores vector frames). Output `Atlas_PDE_Movie.mp4`.

### Shared helper — `reconstruct_stitched_sheet(A, branches, sheet_name, i, j, x_full)`
Single source of truth for rebuilding a sheet's full-domain `[u;q;v;w;p]`
profile at cell `(i,j)`: left half from branch `L` at `(i,j)`, right half from
branch `R` at the mirror cell, spatially flipped. Used by Stages 8, 9, 10 (and
`debug_atlas_point`) so every "grid cell → spatial profile" path is identical.

### Pipeline-local helpers
`file_id` (read a file's `computation_id`), `decide_trace` / `decide_class` (the
§5 staleness checks), `print_provenance` / `print_sheets_provenance` (console
provenance), `run_stage12_movie` (Stage 12 wrapper).

---

## 7. Stability (`p.skip_stability`)

`p.skip_stability = true` (default) skips Phase 5 of the HPSS sweep, leaving
`is_stable_grid` all-false. Setting it `false` runs a per-cell discretised-
Jacobian eigenvalue check on the HPSS branch **only** (HSS and MPSS never compute
stability). It changes nothing unless the HPSS branch is actually rebuilt
(`force_hpss = true`); to propagate the result into `sheet_stable_provisional`
also set `force_sheets = true`. Stability is annotation only and never gates any
existence mask (Invariant 3).

---

## 8. Output files

**Two folder levels.** Every run writes into a per-parameter-set folder
`data/<params>/` where `<params>` encodes the model parameter block
(e.g. `a1_b0p1_c50_d1_e1_f1_Dm0p01_Dc10_L14_MBG28_MP6_res50` — decimal points
become `p`; built by the `atlas_run_name(p)` helper). Within it:

- **`data/<params>/outputs/` + `data/<params>/plots/`** hold outputs that are a
  property of the **parameter set** — the branch libraries, the static atlas,
  and the HSS dispersion stability atlas (Stages 1–7, 2b). Re-running the same
  parameters reuses these.
- **`data/<params>/<sheet>_<i>_<j>/`** is a per-**initial-condition** subfolder
  (e.g. `HSS_HSS_25_26`, named from `trace_sheet` + `trace_cell`) holding every
  **Stage 8–12** output — everything that depends on the PDE trajectory's seed
  sheet + cell. Changing `trace_sheet`/`trace_cell` mints a fresh IC folder, so
  different initial conditions at the **same** parameters never overwrite each
  other while the expensive shared branch libraries/atlas are reused as-is. Both
  folder names are printed at the top of the run.

Stages receive their target via `p.out_dir` (sweeps), or an `out_dir` /
`save_dir` / `save_path` option (every other stage) — param-level stages get the
param folder, Stages 8–12 get the IC folder.

**Param-level data** (`data/<params>/outputs/`, timestamped copies +
`*_latest.mat` canonical): `bvp_sweep_full_results.mat`,
`hss_sweep_multibranch_results.mat`, `branch_{HSS,HPSS,MPSS}_latest.mat`,
`hss_dispersion_atlas.mat`, `atlas_sheets_latest.mat`, plus checkpoints.

**Param-level plots** (`data/<params>/plots/`): branch surface PDFs,
`HSS_Dispersion_Stability.pdf`, `Atlas_Sheet_*.pdf`,
`Atlas_Sheets_Overview.pdf`, `Atlas_Nullclines_Combined.pdf`,
`Atlas_Surfaces_Combined.pdf`.

**Per-IC folder** (`data/<params>/<sheet>_<i>_<j>/`, flat): the trajectory .mat
`atlas_pde_trace_latest.mat` and its classification
`pde_atlas_classification_latest.mat` (+ timestamped copies), the trace-overlaid
phase-plane pages `Atlas_PhasePlane_*.pdf`, the classification summaries
`PDE_Atlas_Classification.pdf` / `PDE_Atlas_Frame_Overlays.pdf`, the movie
`Atlas_PDE_Movie.mp4`, and the `Atlas_PDE_Trace_<id>.png` diagnostic.

---

## 9. Related tools (not part of the automatic sequence)

These tools are **not included in this package**; they live in the development repository.

- `debug_atlas_point(del_pak, del_gef[, sheet])` — inspect existence, stability,
  and the reconstructed profile of every sheet at one grid cell.
- `infer_atlas_transitions` / `plot_atlas_transitions` — segment a classification
  into sustained-dwell transition events.
- `run_atlas_pde_survey` — batch many PDE trajectories from random atlas cells.

---

## 10. Load-bearing rules

- Keep the HSS / HPSS / MPSS distinction explicit; never conflate them.
- MPSS is **constructed** (reflected quarter-domain HPSS), never a direct
  interior-peak BVP.
- Failed cells stay failed (empty + masked + NaN); never back-filled with another
  branch type.
- Stability is annotation only; it never gates existence.
- Classification uses the spatial-profile residual, not phase-plane distance.

