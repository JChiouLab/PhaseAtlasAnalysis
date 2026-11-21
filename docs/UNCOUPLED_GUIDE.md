# Half-cell pipeline: `run_atlas_surface_pipeline.m`

**Role in the paper.** This pipeline computes the steady states of one isolated half cell (Eq. 8 on 0 < x < L/2
with no flux at the middle as at the pole) of the three shapes, over a grid of half-cell amounts. These half-cell
libraries (Stages 1–3) are the starting profiles of the phase-atlas pipeline (`run_seeded_atlas.m`, see
[USER_GUIDE](USER_GUIDE.md)), whose wrapper builds them through `atlas_build_libraries`.

The script's other stages (stitching two isolated halves into sheets, a flux estimated as −2·Dc·(trough difference),
phase-plane plots, simulations classified against the stitched sheets) are an earlier approximation of the atlas.
They are **not used for any figure of the paper**, which computes the coupled quasi-steady states directly.

For the models and their parameters see [MODELS](MODELS.md). The script sits at the repository root; its outputs go
to `data/` (not tracked by git). For the M. xanthus model (`model = 'guzzo_noxa'`) the script builds only the
libraries (Stages 1–3); the parameters are written in the call to `gnx_model('parameters', ...)` near the top of the
script (shipped with `dba_factor` 1; edit the values there).

---

## 1. What the pipeline produces

From one parameter set of the S. pombe model it builds half-cell libraries on the grid of half-cell amounts
`(Δn_PAK, Δn_GEF)` (PAK first, unlike the paper's (ΔGEF, ΔPAK)), stitches them into sheets, then runs one full-model
simulation and compares it with the sheets.

Each half-cell steady state has one of three shapes:

| Library family | Paper | How it is obtained |
|------|---------|--------------------|
| **HSS**  | no-peak (uniform) | algebraic multi-root solve |
| **HPSS** | half-peak (peak at the pole) | a 14-variable BVP on `[0, L/2]` (five concentrations and their derivatives, two sources, two amount accumulators) |
| **MPSS** | full-peak (peak in the interior of the half) | **constructed**: a half-peak solved on a quarter cell, reflected about its peak; not solved directly as an interior-peak problem |

A **sheet** `<L>_<R>` joins a left half of family `L` to a right half of family `R` (the mirror cell). The nine
possible pairings are HSS_HSS, HSS_HPSS, HSS_MPSS, HPSS_HSS, HPSS_HPSS, HPSS_MPSS, MPSS_HSS, MPSS_HPSS and
MPSS_MPSS; which of them exist depends on the parameters. Each sheet carries its own approximate flux field and
nullclines on the shared grid.

**Cell length.** The script ships with `p.L = 7.94`; the paper uses L = 7 (Table 2). Set `p.L = 7` at the top of the
script to match the paper (`M_BG` and `M_P` follow from `p.L`).

---

## 2. Standard usage

```matlab
% edit the parameter and settings block at the top of the file, then:
run_atlas_surface_pipeline
```

A run executes Stages 1–12 in order; each stage either rebuilds its output or reuses the existing file (§5). Every
run writes into its own folder per parameter set, `data/<params>/` (`outputs/` and `plots/`); the outputs of
Stages 8–12, which depend on the simulation's starting point, go into a subfolder `data/<params>/<sheet>_<i>_<j>/`
(§8). Figures are written whether or not windows are shown.

**Typical workflows**

- *First run / new parameter set:* stages whose output file is missing are built anyway. The shipped script has the
  four `force_*` flags `false` (reuse) and `force_trace = true`; set the `force_*` flags `true` to rebuild existing
  libraries and sheets from scratch.
- *Another simulation on an existing set:* leave `force_*` `false` and edit the `trace_*` settings; Stages 8–9 detect
  the change and rerun only the simulation and its classification (§5).
- *Restyle plots only:* edit `trace_n_top_sheets`, `movie_*` or `show_figures`; these never trigger a solve.

---

## 3. Stage sequence

| # | Stage | Function | Output |
|---|-------|----------|--------|
| 1 | half-peak library | `run_hybrid_sweeper_pdefb` → `build_branch_library` | `bvp_sweep_full_results.mat`, `branch_HPSS_latest.mat` |
| 2 | no-peak library | `sweep_hss_multibranch_new` → `build_branch_library` | `hss_sweep_multibranch_results.mat`, `branch_HSS_latest.mat` |
| 2b | linear stability of the uniform state | `plot_hss_dispersion_atlas` | `HSS_Dispersion_Stability.pdf`, `hss_dispersion_atlas.mat` |
| 3 | full-peak library | `build_mpss_branch` | `branch_MPSS_latest.mat` |
| 4 | library surfaces | `visualize_atlas_surfaces` | surface PDFs |
| 5 | stitch sheets | `build_atlas_sheets` | `atlas_sheets_latest.mat` |
| 6 | per-sheet flux plots | `plot_atlas_sheets` | `Atlas_Sheet_*.pdf`, overview, combined nullclines |
| 7 | combined 3-D surfaces | `plot_atlas_surfaces_combined` | `Atlas_Surfaces_Combined.pdf` |
| 8 | full-model simulation | `run_atlas_pde_trace` | `atlas_pde_trace_latest.mat` |
| 9 | classify the simulation | `classify_pde_to_atlas_sheets` | `pde_atlas_classification_latest.mat` |
| 10 | phase-plane pages | `plot_atlas_phase_plane` | `Atlas_PhasePlane_*.pdf` |
| 11 | classification summary | `plot_pde_to_atlas_classification` | `PDE_Atlas_Classification.pdf`, `PDE_Atlas_Frame_Overlays.pdf` |
| 12 | simulation movie | `plot_atlas_pde_movie` | `Atlas_PDE_Movie.mp4` |

Stages 1–3 build the libraries (the part used by the paper); 4–7 display the stitched sheets; 8–12 run and analyse
one simulation against them. Stage 2b (display only) maps the linear stability of the uniform state (§6).

---

## 4. Settings reference (top of the script)

**Model parameters**: `p.a … p.f`, `p.Dm`, `p.Dc`, `p.L`, `p.M_BG`, `p.M_P`, `p.res` (rates, diffusion, cell
length, total GEF and PAK, grid side; e.g. `res = 50` gives 2,500 cells). The mapping `a … f` ↔ Table 2 is in
[MODELS](MODELS.md).

**Solver options** (used only when a library is rebuilt):
- `p.seeding_mode`: `'auto'` (no interaction) or `'interactive'` (pick the first solution visually).
- `p.skip_stability`: `true` skips the per-cell eigenvalue check in the half-peak sweep (§7).
- `p.mpss_strategy`: full-peak construction strategy (`'continuation'`).
- `p.phase0_construct` (default `true`): the half-peak sweep first solves the BVP directly from a constructed edge-peak
  profile, without relaxing a simulation. This reaches a half-peak state that is unstable in time, which relaxation
  cannot reach. If no constructed profile passes the checks, the sweep falls back to relaxation.
- `p.osc_fallback` (default `false` in the solver; the shipped script sets `true`): when the relaxation fallback of
  the half-peak sweep or of the quarter-cell continuation does not settle because the model oscillates, it offers
  `bvp4c` starting profiles from the time mean or from snapshots instead of an arbitrary frame. Needed for parameter
  sets with sustained oscillation. The key setting is `p.pde_T_chunk`: one chunk must span several oscillation
  periods; measure the period with a short `pdepe` run before setting it for new parameters.

**Rebuild flags**: `force_hpss`, `force_hss`, `force_mpss`, `force_sheets` force Stages 1, 2, 3 and 5 to rebuild even
if the output exists (shipped `false`: reuse). `force_trace` forces Stages 8–9 (shipped `true`); normally not needed,
since they detect changes themselves (§5).

**Simulation settings** (Stage 8; a change triggers a new simulation):
- `trace_sheet`, `trace_cell`: the sheet and grid cell `[i_GEF, j_PAK]` the starting profile is taken from.
- `trace_T_end`, `trace_t_fine_end`, `trace_dt_fine`, `trace_dt_coarse`: output times, step `trace_dt_fine` up to
  `trace_t_fine_end`, then `trace_dt_coarse` up to `trace_T_end`.
- `trace_N_mesh`: mesh points over `[0, L]`.

**Classification thresholds** (Stage 9): `trace_ambiguity_gap`, `trace_off_atlas_tol` flag a frame as ambiguous or
off the sheets. They are not calibrated; the residuals themselves are the main output.

**Display only** (never trigger a solve):
- `trace_n_top_sheets` (Stage 10): `[]` draws one page per sheet the simulation stayed on; `N` draws the N
  best-fitting sheets; `Inf` or 9 draws all nine.
- `render_movie`, `movie_fps`, `movie_frame_stride`, `movie_resolution_scale`, `movie_quality` (Stage 12).
- `show_figures`: `true` shows figure windows; `false` runs without windows. Every figure is exported to
  `data/<params>/plots/` either way.

---

## 5. Rebuild and reuse

**Stages 1–3 and 5** rebuild only when their `force_*` flag is `true` or their output `.mat` is missing. Rebuilding a
library changes its `computation_id`.

**Stages 8–9** follow the settings. Stage 8 stores the settings that define the simulation (`trace_sheet`,
`trace_cell`, the output times, `trace_N_mesh` and the `computation_id`s of the libraries and sheets it used) in the
trace file as `run_config`; on the next run `decide_trace` compares them and reruns on any difference (including a
rebuilt library). Stage 9 (`decide_class`) reclassifies when the simulation was rerun, the classification is missing
or was made from another simulation, or the thresholds changed.

**Stages 4, 6, 7, 10, 11, 12** always run; they only display computed data.

---

## 6. Function reference

### Stage 1: `run_hybrid_sweeper_pdefb(p)`
Solves the half-peak BVP at every grid cell at fixed `L`. The first solution, at the grid centre, comes from the
constructed edge-peak profile (`p.phase0_construct`) or else from relaxation of a simulation (with a gradual change of
L if needed). It is then continued across the grid by predictor–corrector; where that fails, a neighbour's profile is
relaxed in a simulation (oscillation-aware if `p.osc_fallback`) and given to `bvp4c`. Failed cells are masked, never
filled in. An optional last step (unless `p.skip_stability`) computes each solved cell's stability in time from the
eigenvalues of the discretized Jacobian.

Every converged solution must pass the shape classifier `classify_hpss_morphology` before it is stored or used as a
starting profile. A valid half-peak may have a slight dip at the pole and its maximum within the first 10 % of the
half, but the other end must stay the minimum, with no second significant maximum or interior minimum. The stored
solutions are checked again before the recovery passes and after each of them; solutions that collapsed to the
uniform state, interior-peak, multi-peak or unclassified profiles are set aside as failed half-peak cells. Recovery
can therefore approach the edge of the half-peak region but cannot carry another shape across it.

Output `bvp_sweep_full_results.mat`: `PAK_DIFF`, `GEF_DIFF`, `{U,Q,V,W,P}_trough`, `SOL_grid` (BVP solutions),
`failure_mask`, `is_stable_grid`, geometry, `par`, `p`, `computation_id`, `morphology_class_grid`,
`morphology_reason_grid`, `morphology_diagnostics_grid` and `quarantined_nonhpss_mask` (cells where a converged
solution of another shape was set aside; it does not mean that no half-peak state exists there).

### Stage 2: `sweep_hss_multibranch_new(p)`
Finds all real positive uniform steady states per cell by `fsolve` from many starts, then labels the branches by
continuity in `[u q v w p]` between neighbouring cells (not by the order of the roots in a cell), so one branch keeps
one label across the grid. Output `hss_sweep_multibranch_results.mat`: `{U,Q,V,W,P}_hss_branches` (`res×res×maxK`),
`n_branches_grid`, `branch_id_grid`, `branch_summary`.

### Stages 1–2 adapter: `build_branch_library(hpss_path, hss_path, out_dir)`
Reads the half-peak and no-peak sweep files and writes them in the common library format used downstream. A cell
counts as half-peak only if the solver converged **and** the profile is a half-peak at the pole (amplitude ≥ 0.5),
by the same classifier as the sweep. Writes `branch_HPSS_latest.mat` and `branch_HSS_latest.mat` with per-cell
`branch_class`, `morphology_reason`, `existence_mask`, `SOL`, `flux_values`, `flux_source_location`, `amp`, `peak_x`,
`mass_residual` and stability fields (information only); never modifies the source files.

### Stage 2b: `plot_hss_dispersion_atlas(opts)`
Linear stability of the uniform state in each cell from the dispersion relation: the reaction Jacobian `J` (analytic;
it satisfies both conservation laws, null vectors `[1 1 1 1 0]` and `[0 1 0 0 1]`) is scanned over the no-flux modes
`kₙ = nπ/ℓ` (`ℓ = L/2`, `n ≥ 1`), and the cell is marked unstable when `max Re eig(J − kₙ²D) > tol`, with the
fastest-growing mode `n*`. Two panels over `(Δn_PAK, Δn_GEF)`: (A) the unstable region, coloured by `n*/2` periods on
the half cell (0.5 = one half-peak); (B) whether the onset is stationary or oscillatory. `n*/2` is the linear onset,
not the number of peaks of the final state (peaks of mass-conserving systems coarsen). Options: `hss_path`,
`save_dir` (required), `domain` (`'half'` default or `'full'`), `n_max` (80), `tol` (1e-6), `visible`. Writes
`HSS_Dispersion_Stability.pdf` and `hss_dispersion_atlas.mat`.

### Stage 3: `build_mpss_branch(p)`
Builds the full-peak library: a half-peak sweep on a quarter cell, then reflection about the peak into a half-cell
interior-peak profile. The flux is taken at the outer minimum, at the middle of the cell. Each grid row and each
rescue cell runs with a wall-clock limit (`p.phase3_row_timeout`, `p.cell_timeout`) and checkpoints, so the build can
resume; a row or cell that times out is left failed for the rescue pass. Failed quarter-cell solves or reflections
are masked. Output `branch_MPSS_latest.mat` in the library format (stability is not computed for this shape).

### Stage 4: `visualize_atlas_surfaces(opts)`
Plots the three library files: existence masks, 3-D flux surfaces (`v+w` and `p`, mirrored through the origin) and
scalar fields (`amp`, `peak_x`, `mass_residual`). Nullclines are added only if a sheets file exists. Display only.

### Stage 5: `build_atlas_sheets(hpss_path, hss_path, mpss_path, out_dir)`
Stitches the nine pairings from the three library files, with no new solve. The mirror of cell `(i,j)` is exactly
`(res+1-i, res+1-j)` (index arithmetic on the symmetric grid, no interpolation). Per sheet:
`sheet_exists = left.existence(i,j) & right.existence(mirror)`; approximate fluxes
`U_flux = -2·Dc·(p_left − p_right)` and `V_flux = -2·Dc·(vw_left − vw_right)` from each library's stored
`flux_values`; `sheet_stable_provisional` (information only). This flux is not the paper's Eq. (9). Output
`atlas_sheets_latest.mat`: a `sheets` struct of the nine pairings plus `sheet_names`, geometry, `par` and the
`computation_id`s of the sources.

### Stage 6: `plot_atlas_sheets(opts)`
One page per sheet: existence region, flux arrows of equal length and the two nullclines, from Stage 5's stored
flux. Also a 3×3 overview and a combined nullcline figure.

### Stage 7: `plot_atlas_surfaces_combined(opts)`
One 3-D figure of the three libraries' flux surfaces, with every sheet's nullclines drawn at their height (on a
nullcline the left and mirrored right surfaces have equal height, so the height is unique).

### Stage 8: `run_atlas_pde_trace(opts)`
Builds a whole-cell starting profile from the chosen sheet (with `reconstruct_stitched_sheet`, so the first frame
equals the classifier's reconstruction), integrates the full model with `pdepe` over `t_span`, and records the
imbalances `(Δn_PAK(t), Δn_GEF(t))`. Output `atlas_pde_trace_latest.mat`: `sol_pde` (`nt×nx×5`), `x_full`,
`t_span`, `nPAK_trace`, `nGEF_trace`, `u0_stitched`, the start and `run_config`.

### Stage 9: `classify_pde_to_atlas_sheets(opts)`
For every frame, takes the nearest grid cell and computes, for each sheet that exists there, the relative L2
difference between the frame's profile and the sheet's reconstructed profile (a profile difference, not a distance
in the plane). The best sheet has the smallest difference; frames are flagged ambiguous (second best within
`ambiguity_gap`) or off the sheets (best above `off_atlas_tol`). Output `pde_atlas_classification_latest.mat`:
`t_span`, `dpak`, `dgef`, `best_cell`, `residual_by_sheet` (`nt×9`), `residual_by_species`, `best_sheet_idx/name`,
`best/second_best_residual`, `ambiguous`, `off_atlas`.

### Stage 10: `plot_atlas_phase_plane(opts)`
One page per selected sheet with its flux field and nullclines, the whole simulation as a grey dotted line and the
parts on this sheet as a thick black line. Needs Stage 9. Output `Atlas_PhasePlane_<sheet>.pdf`.

### Stage 11: `plot_pde_to_atlas_classification(opts)`
The best sheet per frame over time, the profile difference to all nine sheets over time, and an overlay of frames
and their best sheet profiles (all five species at a few times). Outputs `PDE_Atlas_Classification.pdf` and
`PDE_Atlas_Frame_Overlays.pdf`.

### Stage 12: `plot_atlas_pde_movie(opts)`
Renders the simulation to `.mp4`: two panels (`u, q` and `v, w, p` against `x`), one video frame per stored frame
(every `frame_stride`-th), fixed colours and axis limits; `resolution_scale` supersamples the raster. Output
`Atlas_PDE_Movie.mp4`.

### Shared helper: `reconstruct_stitched_sheet(A, branches, sheet_name, i, j, x_full)`
Builds a sheet's whole-cell `[u;q;v;w;p]` profile at cell `(i,j)`: the left half from family `L` at `(i,j)`, the
right half from family `R` at the mirror cell, flipped. Used by Stages 8, 9 and 10.

### Helpers in the script
`file_id` (reads a file's `computation_id`), `decide_trace` / `decide_class` (the checks of §5),
`print_provenance` / `print_sheets_provenance` (console output), `run_stage12_movie` (Stage 12).

---

## 7. Stability (`p.skip_stability`)

`p.skip_stability = true` (default) skips the stability step of the half-peak sweep, leaving `is_stable_grid`
all false. With `false`, each half-peak cell gets an eigenvalue check of the discretized Jacobian (no-peak and
full-peak cells never do). It has an effect only when the half-peak library is rebuilt (`force_hpss = true`); to carry
it into `sheet_stable_provisional`, also set `force_sheets = true`. Stability is information only; it never decides
whether a cell exists.

---

## 8. Output files

**Two folder levels.** Every run writes into a folder per parameter set, `data/<params>/`, where `<params>` lists the
parameters (decimal points become `p`; built by `atlas_run_name(p)`), e.g.
`a0p595_b0p185_c75p8_d0p552_e0p828_f1_Dm0p05_Dc5_L7_MBG14_MP3p444_res50`. Inside it:

- **`outputs/` and `plots/`** hold what depends only on the parameters: the libraries, the stitched sheets and the
  stability map (Stages 1–7, 2b). Running the same parameters again reuses them.
- **`<sheet>_<i>_<j>/`** (e.g. `HSS_HSS_25_26`, from `trace_sheet` and `trace_cell`) holds every output of
  Stages 8–12. Another starting point gets its own folder, so simulations at the same parameters never overwrite each
  other. Both folder names are printed at the start of the run.

Data (`outputs/`, with time-stamped copies and `*_latest.mat`): `bvp_sweep_full_results.mat`,
`hss_sweep_multibranch_results.mat`, `branch_{HSS,HPSS,MPSS}_latest.mat`, `hss_dispersion_atlas.mat`,
`atlas_sheets_latest.mat`, plus checkpoints.

Plots (`plots/`): library surface PDFs, `HSS_Dispersion_Stability.pdf`, `Atlas_Sheet_*.pdf`,
`Atlas_Sheets_Overview.pdf`, `Atlas_Nullclines_Combined.pdf`, `Atlas_Surfaces_Combined.pdf`.

Per starting point (`<sheet>_<i>_<j>/`): `atlas_pde_trace_latest.mat`, `pde_atlas_classification_latest.mat`
(with time-stamped copies), `Atlas_PhasePlane_*.pdf`, `PDE_Atlas_Classification.pdf`, `PDE_Atlas_Frame_Overlays.pdf`,
`Atlas_PDE_Movie.mp4` and `Atlas_PDE_Trace_<id>.png`.

---

## 9. Rules the code keeps

- No-peak, half-peak and full-peak libraries are kept apart; a failed cell stays failed (empty, masked, NaN) and is
  never filled with another shape.
- The full-peak half is constructed from a reflected quarter-cell half-peak, never solved directly.
- Stability is information only; it never decides existence.
- Simulations are compared with the sheets by profile difference, not by distance in the plane.
