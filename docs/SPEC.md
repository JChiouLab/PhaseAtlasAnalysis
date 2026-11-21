# Design and guarantees

What a run of the phase-atlas pipeline (`run_seeded_atlas.m`) guarantees, and what it does not. Usage is in the
[user guide](USER_GUIDE.md); the models and tested settings are in [MODELS](MODELS.md).

## Design

- **Two pipelines.** The half-cell pipeline (`run_atlas_surface_pipeline.m`) builds libraries of isolated half-cell
  steady states. The phase-atlas pipeline reads them as starting profiles and solves the whole-cell problem of
  Eq. (8) at every grid point of the plane. Each model sits behind a model module (equations, fields, conserved
  quantities, scaling, library reading, model-specific checks); the shared code names no model constant.
- **Routes to a solution.** Cold starts from two library halves, continuation from neighbouring solutions, a per-root
  pass to every neighbour, rescue of unreached points, and, optionally, solutions of an earlier run (`import`;
  `focus` mode uses only chosen ones and no library). All routes write to one store of solutions.
- **Modes.** `smoke`, `pilot` and `production` differ in where cold starts are placed, how far continuation goes and
  the default budget; `focus` refines chosen sheets of an earlier run (USER_GUIDE section 4).
- **Restart.** Every call is written to its own file as soon as it returns. A stopped or killed run resumes by
  replaying the stored calls without solving and continuing; the effective settings are saved with the run.

## What a run guarantees

- **Every stored solution passed the same checks**, at tolerances fixed before the run: finite values, the residual
  of the equations, the boundary conditions and continuity at the middle, the amounts in each half, positivity
  everywhere, and agreement between the flux given by the sources and the flux read from the profile slopes at the
  middle. A small solver residual alone is not enough. Calls that fail are kept with the reason and never used as
  starting points.
- **One solution per distinct profile.** Two solutions at the same point count as one when their profiles differ by
  less than 1e-3 (scale-aware distance over all fields); distinct solutions are all kept. The family label names the
  shape and never creates or removes a solution.
- **Traceable settings.** `settings_effective.txt` lists every setting with its value and source; unknown or
  contradictory settings are refused before any call.
- **Bounded cost.** A run stops at its call budget or wall time and can be resumed; refused, unattempted and
  policy-skipped points are reported separately.
- **Same code, same records.** A run records the software version, the model code, the MATLAB release (with the
  Update), the architecture, the thread mode and the solver file hashes. A resume refuses any difference instead of
  mixing them.
- **Import is never a hidden start.** Imported solutions are re-solved and checked like any other call; the source
  run is identified by a hash and recorded in the summary.

## What a run does not guarantee

- **Completeness.** A point without a solution is not proof that none exists; only the configured search was done.
  Under the default production targets, coexisting solutions at a point are looked for only by continuation and the
  per-root pass (the summary says so).
- **Stability.** No stored solution is shown to be stable in time.
- **Branch connectivity.** Neighbouring solutions are not shown to lie on one continuous branch.
- **Simulations.** A full-model simulation (trace) is an integration from one solution; its end state is not shown to
  be an attractor. Plots are displays, not verified results.
- **Other settings or releases.** The tolerances, mesh sizes and shape thresholds were chosen at the tested settings
  (MODELS). Another release of MATLAB is accepted after a self-check of the patched `bvp4c`, which tests the patch,
  not the numerics: it may find other solutions, call counts or meshes.

## Long cells and higher modes

- `seed_modes` chooses the half shapes used for cold starts (default: the paper's three); more shapes come from
  higher-mode libraries built by `atlas_build_libraries(..., struct('max_mode', n))`. An empty mode library means the
  search found nothing, not that the mode does not exist.
- Above a reference length (S. pombe 14 µm, M. xanthus 6 µm) the mesh and node counts scale with L (USER_GUIDE
  section 6). At L = 30 the scaled cold mesh cap cut most higher-mode cold calls; a root not found under the cap is a
  search result, and `cold.mesh_cap` is the setting to change.
