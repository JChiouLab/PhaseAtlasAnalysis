# PhaseAtlasAnalysis

The code of the paper "Phase atlas geometry of reaction-diffusion systems underlies cell polarity transitions"
(Wang, Saulnier, Hinds, Hsu, Lin, Chiou). It computes the **phase atlas** of a mass-conserving polarity model: the
cell is cut into two halves at its middle, and at each point of the plane of mass imbalances between the halves the
pipeline finds the coupled quasi-steady states (QSS) of the two halves (Methods, "Phase-plane construction").

Two MATLAB pipelines and one Python package:

- **the phase-atlas pipeline** (`run_seeded_atlas.m`): solves the whole-cell boundary value problem of Eq. (8) at
  every grid point of the plane, starting from half-cell profiles and then from neighbouring solutions, and stores
  every quasi-steady state with its sources; it also plots the states and runs full-model simulations started from
  them;
- **the half-cell pipeline** (`run_atlas_surface_pipeline.m`): computes the steady states of one isolated half cell
  (no flux at the middle) of the three shapes over a range of amounts; these are the starting profiles of the
  phase-atlas pipeline (see `docs/UNCOUPLED_GUIDE.md`);
- **`noisy_simulations/`** (Python): the noisy simulations of Fig. 4 (see its `README.md`).

Version: see `VERSION` and `CHANGELOG.md`. Licence: CC BY-NC 4.0 (non-commercial use with attribution; see
`LICENSE`).

## Paper figures

Which computation underlies each figure. The code's names for the paper's quantities are in the appendix
("Paper ↔ code").

- **Fig. 2, 2E, S2 (M. xanthus):** the planes in the paper were computed with the authors' first M. xanthus
  implementation of the same coupled boundary value problem (its 24-state form is kept here as
  `models/guzzo_noxa/sag_coupled24.m`). This pipeline computes the same planes (Example 1 below) and was checked
  against them. Δn2 of
  Fig. S2D is computed on membrane MglB (`Bm`, the default).
- **Fig. 3, S3, S4, Movie S1 (S. pombe, L = 7):** the phase-atlas pipeline with `model = 'original_polarity'`, the
  parameters shipped in `run_seeded_atlas.m` (L = 7, MG = 14, MP = 2 · 0.246 · 7 = 3.444), `run_mode =
  'production'` and a 101 × 101 grid (Example 2 below). The limit cycles are full-model simulations
  (`seeded_atlas_trace`) started from quasi-steady states of that run.
- **Fig. 5, S5 (Ashbya setting, L = 30):** the same pipeline and model with L = 30, `M_BG = 60` and
  `M_P = 2*(0.246/3)*30` (= 4.92), a 101 × 101 grid, as `focus` runs on chosen sheets seeded from an earlier
  51 × 51 run (Example 3 below).
- **Fig. 4:** `noisy_simulations/` (Python).

Full-model simulations started from a quasi-steady state are run by `seeded_atlas_trace` (ode15s for S. pombe, pdepe
for M. xanthus). Other simulation set-ups are not included.

## Requirements

- MATLAB with the Java VM (do not start MATLAB with `-nojvm`).
- Tested on MATLAB R2025b Update 4, macOS arm64. Other releases are accepted after a self-check (below) but are
  untested; the oldest release that works is not known (the functions used date from R2020a or earlier).
- Toolboxes: Optimization (`fsolve`), Signal Processing (`findpeaks`) and Parallel Computing (`parfor`, `parfeval`)
  are needed only to build the S. pombe half-cell library and to run the half-cell pipeline's other stages. The
  M. xanthus library, the phase-atlas runs, plots and simulations need base MATLAB only.
- The 2-D sheet PDFs: Python 3 with numpy, matplotlib and contourpy
  (`src/pipelines/seeded_atlas_standalone/post/requirements_sheets.txt`). MATLAB runs `python3` from the system path;
  set the environment variable `SA_PYTHON` to use another. Without them the runs and all other outputs are unchanged,
  and the command that draws the PDFs later is printed.
- `noisy_simulations/`: Python 3 (see its `requirements.txt`).

## Setup

From the repository folder in MATLAB:

```matlab
setup_atlas            % puts src/pipelines/seeded_atlas_standalone on the path and makes the traced bvp4c copy
setup_atlas('check')   % verifies only; prints one line, OK or FAIL with the reason
```

The pipeline records the mesh sizes `bvp4c` proposes. For this, `setup_atlas` copies **your own MATLAB's**
`bvp4c.m` into `src/pipelines/seeded_atlas_standalone/solver/sa_bvp4c_traced.m` (and `solver/private/`) with a
renamed function line and two added lines, then runs a self-check of a few seconds (about 12 small `bvp4c` calls).
These files are not tracked by git; no MathWorks file is in this repository. `run_seeded_atlas.m` calls
`setup_atlas('quiet')` itself before every run.

## Quick start

1. **A small run.** Open `run_seeded_atlas.m`, edit the block at the top and run the file. For a first try set
   `run_mode = 'smoke';` and a small library grid, e.g. S. pombe `res = 15;` (or the M. xanthus example in the
   block's comment with `res = 11;` instead of 25). When the library folder is empty the script first builds the
   half-cell library (the S. pombe build opens a parallel pool and can take an hour or more). It prints every path
   before any work.
2. **Plot.** The run saves `plots/roots_3d.fig` (the atlas in three dimensions, rotatable) and the 2-D sheets as
   PDFs in the style of Fig. 5, `plots/sheets_2d/<sheet>.pdf` and `all_families.pdf` (existence, flux arrows,
   nullclines and an example profile); redraw with `seeded_atlas_plot(run_dir)`. *Note: Root family identification is based on arbitrary
  thresholds and can be wrong. Manual inspection of the peak locations my be
  needed in some cases.*
3. **One simulation.** A full-model simulation started from one quasi-steady state:
   `seeded_atlas_trace(run_dir, struct('name', 'T1', 'root', '<a root_id from roots_index.csv>'))`, or list it in
   the script's `traces` block.
4. **A half-cell library on its own.** `setup_atlas` adds only the phase-atlas package, so first add the half-cell
   pipeline:
   ```matlab
   addpath('src/solvers/atlas_surface_pipeline_dependencies');
   % M. xanthus also: addpath('src/solvers/atlas_surface_pipeline_dependencies/models/guzzo_noxa'); addpath('src/models/guzzo_noxa');
   par = struct('L', 6, 'r_speed', 0.35, 'dba_factor', 0.35, 'rho_R', 1, 'rho_B', 1, 'dn_half_range', 1.8);
   m = atlas_build_libraries('guzzo_noxa', par, 25, fullfile(pwd, 'data', 'my_library'));
   ```
5. **The half-cell pipeline:** edit and run `run_atlas_surface_pipeline.m` (see `docs/UNCOUPLED_GUIDE.md`).
6. **Starting from an earlier run.** Set `import_from` to a finished run folder of the same model and parameters:
   its quasi-steady states (and their mirror images) are solved once at the nearest points of the new grid and
   serve as starting solutions. `run_mode = 'focus'` uses only the states you pick and restricts the search to their
   sheets. See `docs/USER_GUIDE.md` section 3.

Documents: `docs/USER_GUIDE.md` (the phase-atlas pipeline: settings, modes, results, resume),
`docs/UNCOUPLED_GUIDE.md` (the half-cell pipeline), `docs/MODELS.md` (the two models, their parameters and the tested
settings), `docs/SPEC.md` (what a run guarantees and what it does not).

## Examples: the paper's three atlases

Three set-ups, from simple to involved. Each is written as you would run it with this version of the code. Start
MATLAB in the repository folder and run `setup_atlas` once. Each run prints its run folder before any work, and
running the same set-up again resumes a stopped run.

### Example 1: a simple atlas (M. xanthus, Fig. 2C)

The M. xanthus plane at signal factor f = 0.02. Every setting except the parameters and the call limit is the
default: production mode, the paper's three half shapes, and the model's grid step of 0.15 in `dn`
(27 × 27 = 729 grid points). The call is the one `run_seeded_atlas.m` makes after its edit block.

```matlab
setup_atlas;
cfg = struct('model', 'guzzo_noxa', 'mode', 'production', 'res', 25, ...
    'par', struct('L', 6, 'r_speed', 0.35, 'dba_factor', 0.02, 'rho_R', 1, 'rho_B', 1, 'dn_half_range', 1.8), ...
    'budget', struct('max_calls', 30000, 'max_wall_s', 36000));
W = seeded_atlas_wrapper(cfg);                          % builds the half-cell library if missing, then solves
```

The run takes about 1,700 solver calls and 8 minutes and finds 431 quasi-steady states. The half-cell library build
needs base MATLAB only. The results are in
`data/standalone_runs/guzzo_noxa_production_L6_rspeed0p35_dbafactor0p02_rhoR1_rhoB1_dnhalfrange1p8_res25_step0p15_<hash>/`:
`roots_index.csv` lists every state with its coordinates and sheet, `plots/roots_3d.fig` is the 3-D atlas and
`plots/sheets_2d/` holds the sheet PDFs. For the other signal levels of Fig. 2E, change `dba_factor`.

### Example 2: a large atlas with simulations (S. pombe, L = 7, Fig. 3, S3, S4)

The Fig. 3 atlas on a 101 × 101 grid, with the two simulations drawn in Fig. 3. In `run_seeded_atlas.m` keep the
shipped `model` and `par` (L = 7, MG = 14, MP = 3.444) and set:

```matlab
res      = 50;                                  % half-cell library grid (built on first use; opens a parallel pool)
points_per_axis = 101;                          % 101 x 101 grid points over the plane
run_mode = 'production';
budget   = struct('max_calls', 60000, 'max_wall_s', 7*86400);
settings_overrides = struct();
settings_overrides.targets.policy = 'all_seedable';   % a cold start at every grid point a library pair covers
traces(1).name = 'LC1';  traces(1).root = 'k:-2,-17/R000001';   % settles on the limit cycle
traces(2).name = 'LC2';  traces(2).root = 'k:-5,-10/R000001';   % damped oscillation into a polarized state
```

- **`targets.policy = 'all_seedable'`** starts the solver from the half-cell library at every grid point where a
  pair of library halves exists (about 20,000 cold starts here). This finds more of the points where several states
  coexist. The default, `'per_seed'`, makes one cold start per library entry (about 5,200 here). It is faster but
  can miss some of those states; the paper's atlas was made with cold starts at every point.
- **Budget and resume.** The paper's run took 40,527 solver calls in six sittings and found 5,461 states. When a run
  stops at `max_wall_s`, run the file again; after a `max_calls` stop, raise `max_calls` first.
- **Traces.** A trace is a full-model simulation (ode15s, to t = 5000 by default) started from a quasi-steady state.
  The two root ids name grid points of this run (`k:<i>,<j>`, steps from the centre). Choose your own from
  `roots_index.csv` or by clicking on `plots/roots_3d.fig` (`traces(1).at = [x y z]`). Traces listed after a run has
  finished are integrated on the next run of the file, with no new solves.

### Example 3: refining chosen sheets with a focus run (Ashbya setting, L = 30, Fig. 5, S5)

At L = 30 the plane holds many sheets, and a full 101 × 101 production run would spend most of its calls on sheets
the figure does not show. The paper's Fig. 5 sheets were therefore made in two steps.

**Step 1: a coarse atlas of every sheet.** A production run on 51 points per axis. Profiles are longer in a
30 µm cell, so the solver's mesh limits are raised in proportion to L (10,000 × 30/7):

```matlab
par      = struct('a', 2.38/4, 'b', 0.185, 'c', 151.6/2, 'd', 0.552, 'e', 0.828, 'f', 1, ...
                  'Dm', 0.05, 'Dc', 5, 'L', 30, 'M_BG', 2*30, 'M_P', 2*(0.246/3)*30);   % M_P = 4.92
res      = 25;
points_per_axis = 51;
run_mode = 'production';
budget   = struct('max_calls', 50000, 'max_wall_s', 7*86400);
settings_overrides = struct();
settings_overrides.solver.nmax        = 42858;  % mesh points allowed per call
settings_overrides.solver.mesh_reject = 42858;
settings_overrides.cold.mesh_cap      = 42858;  % mesh cap of the cold starts from the library
```

The paper's step-1 run took 29,256 calls (about 3.6 days) and found 7,743 states. Continuation reached many more
sheets than the three starting shapes give, including halves with several peaks.

**Step 2: focus runs on the chosen sheets.** A focus run takes the states of the step-1 run on the sheets you
choose, solves each once at the nearest point of a finer grid, and continues from them. It makes no cold starts, does
not need the half-cell library, and does not step beyond the chosen sheets and their mirror images.

```matlab
run_mode    = 'focus';
import_from = '<the step-1 run folder>';
Z = load(fullfile(import_from, 'internal', 'settings_effective.mat'), 'parameters');
par         = Z.parameters;                     % the source's parameters exactly (import requires them bit for bit)
points_per_axis = 101;                          % twice as fine as the source
budget      = struct('max_calls', 20000, 'max_wall_s', 7*86400);
settings_overrides = struct();
settings_overrides.import.families = {'1+/0', '0/1-'};   % half-peak/no-peak and its mirror
run_dir     = fullfile(pwd, 'data', 'standalone_runs', 'L30_focus_n101_1p0');
```

Fig. 5 used three such runs, one per sheet and its mirror: `{'1+/0','0/1-'}`, `{'2-/0','0/2-'}` and
`{'1+/2-','2-/1-'}`. Fig. S5 also uses `{'0/0','1+/1-'}`. Each focus run solves one call at a time, so the
groups can run at the same time in separate MATLAB processes. For example, save one copy of `run_seeded_atlas.m` per
group with the block above and launch each from a terminal:

```
matlab -singleCompThread -noFigureWindows -batch "run('my_L30_focus_g1.m')" > g1.log 2>&1 &
```

Use the same launch flags when you resume. Section 3 of `docs/USER_GUIDE.md` ("Focus mode") lists the focus
settings: choosing single states with `import.root_ids`, the family fence `continuation.families`, and reuse of
the source's failed steps.

## Folder map

```
setup_atlas.m                    path + traced bvp4c copy + self-check
run_seeded_atlas.m               the phase-atlas pipeline (edit the block at the top)
run_atlas_surface_pipeline.m     the half-cell pipeline
src/pipelines/seeded_atlas_standalone/              the phase-atlas package (models/, post/, solver/, private/)
src/solvers/atlas_surface_pipeline_dependencies/    the half-cell pipeline (incl. models/guzzo_noxa/)
src/models/guzzo_noxa/                               the M. xanthus model equations, shared by both pipelines
tests/seeded_atlas_standalone/   the tests: run_portable_tests
docs/                            USER_GUIDE, UNCOUPLED_GUIDE, MODELS, SPEC
noisy_simulations/               the Fig. 4 noisy simulations (Python)
data/                            created by the first run; not tracked by git
```

**Outputs** go to `data/` in the repository folder (half-cell libraries in `data/<parameters>/outputs`, phase-atlas
runs in `data/standalone_runs/<name>`), unless you give a folder.

**Tests:** `addpath('tests/seeded_atlas_standalone'); run_portable_tests` (about 15 minutes; mostly stand-in models
plus about 140 small `bvp4c` calls). It prints PASS, FAIL or SKIP per test and writes scratch files to
`data/experiments/portable_tests_<stamp>/`.

## Other MATLAB releases

The pipeline accepts any release whose `bvp4c.m` can be copied with the two added lines and whose copy passes the
self-check; otherwise it stops with a message naming the failed check. On a release that passes, the solver is that
release's own `bvp4c`, and the pipeline code, settings and acceptance checks are the same. The same quasi-steady
states, call counts, meshes or timings as on R2025b are not guaranteed, nor are bit-identical results across
releases or machines. A run records the MATLAB release, its Update, the architecture and the toolbox file hashes,
and refuses to resume on another. The half-cell pipeline uses your own `bvp4c`, `pdepe` and `fsolve` without a
self-check.

## Platform notes

1. **Line endings.** `.gitattributes` keeps every file byte-exact on every platform; code hashes depend on it. Do not
   re-save files with an editor that converts line endings.
2. **Windows (untested).** Default run-folder names are about 90 characters, and paths under
   `data/standalone_runs/<name>/` can exceed 260 characters: clone to a short path or enable long paths.
3. **`metadata_log.csv`** appears in the current folder after an S. pombe library build (not tracked by git).
4. **Pool size:** open `parpool(n)` first or pass `opts.pool_workers`; the default pool uses all cores. One worker
   gives no parallelism but runs two MATLAB processes.
5. **Tests during a run:** do not run the tests while a run is in progress in the same repository folder.
6. **MATLAB Updates:** a run cannot be resumed after a MATLAB Update. Finish it first, or start a new `run_dir`.

## Known limitations

- **The half-peak library can hold a few unconverged entries.** The half-peak sweep stores the solver's last
  iterate without a convergence check, so some stored cells have a residual above the solver's tolerance (in the
  authors' L = 7 library, 12 of 1,096). As starting profiles this is harmless, since every coupled solution is
  checked on its own, but existence read from the half-cell masks is overstated.
- **The full-peak and higher-mode constructors check the slope only at the ends of the joined pieces**, where the
  no-flux condition makes it zero, so the slope at the interior join is not checked.
- Windows is untested; other MATLAB releases are self-checked, not tested.
- No resume across a MATLAB Update or another architecture.
- A `bvp4c` that cannot be copied and patched (a `.p` or MEX file, or a rewritten one) is refused.

## Provenance

Exported from the authors' development repository; history not included. `EXPORT_MANIFEST.txt` lists every file
with its SHA-256. The user guide keeps its section numbers, which the code cites as "USER_GUIDE section N".

`src/pipelines/seeded_atlas_standalone/post/guzzo_noxa/Guzzo_noXA.m` is the reference script of the reduced
M. xanthus model, kept byte for byte. It is never executed; its hash is part of the identity of M. xanthus
simulations.

## Credits

Pipeline and documentation: J.-G. Chiou, with AI assistance as declared in the paper. Reduced M. xanthus model and
its reference script `Guzzo_noXA.m`: W.-H. Wang. Implementation of the M. xanthus model in the pipeline: JB
(Jean-Baptiste Saulnier). Noise package `gef_pak_noise`: C. Hinds.

## Appendix: paper ↔ code

| paper | code |
|---|---|
| the M. xanthus model (Eqs. 1-4, Table 1) | `guzzo_noxa` (Guzzo et al. 2018 without MglA: "noXA"); fields `Rm Rc Bm Bc` = rm rc bm bc; `Aeff` = a_eff; `KA` = Ka; `DM`, `DC` = Dm, Dc |
| the S. pombe (fission yeast) model (Eq. 5, Table 2) | `original_polarity`; fields `u q v w p` as in the paper |
| no-peak / half-peak / full-peak half | library families `HSS` / `HPSS` / `MPSS`; shape labels `0` / `1+` (left half) or `1-` (right half) / `2-` |
| sheet, named left/right | family label `<left>/<right>`: `1+/0` = half-peak/no-peak, `0/1-` = no-peak/half-peak, `2-/0` = full-peak/no-peak, `2-/1-` = full-peak/half-peak, `1+/1-` = half-peak/half-peak, `0/0` = no-peak/no-peak; `1-/1+` = one peak at mid-cell |
| quasi-steady state at a point of the plane | a *root* at a lattice *node* (grid point) |
| starting profiles from isolated halves | the *seed library* (`branch_HSS/HPSS/MPSS_latest.mat`) |
| sources s_L, s_R of Eq. (8) | the shifts (`sigma`, `s_GL` ...) stored in each solution; the flux J of Eq. (9) follows from them |
| S. pombe coordinates ΔGEF/MG, ΔPAK/MP | `d_G` = ΔGEF/MG; `d_P` = ΔPAK/MG, so ΔPAK/MP = `d_P`/ρ |
| M. xanthus coordinates ΔRomR/MR, ΔMglB/MB | `dn_R`, `dn_B` = (n_left − n_right)/(L/2), a density difference; Δn/M = `dn`/(2ρ) with ρ = `rho_R` or `rho_B` (1 in the paper) |
| shape axis Δn2 (Eq. 10) | `dn2` = ∫₀¹ cos(2πη) m dη on the rescaled cell η = x/L, in the model's own concentration units; m = `u` (S. pombe) or `Bm` (M. xanthus) by default (`model_options.main_variable`). M. xanthus: Δn2 = L · `dn2`. S. pombe (concentrations scaled so that MG = 1): Δn2 = MG · `dn2` = 2L · `dn2` |
| Table 2: kgg, dg, kpg, dgp, k*, kg | `a`, `b`, `c`, `d`, `e`, `f` (unrounded: `2.38/4`, `0.185`, `151.6/2`, `0.552`, `0.828`, `1`) |
| Table 2: MG = 2L, MP = 2ρL | `M_BG = 2*L`, `M_P = 2*0.246*L`: the code uses ρ = 0.246, which Table 2 rounds to 0.25; Ashbya set ρ = 0.246/3 |
| Table 1: kr, krb, dr | `r_speed` × `kR`, `kRB`, `dR`; Table 1 lists the products at `r_speed = 0.35` |
| signal × f | `dba_factor` (multiplies `dBA`) |
| full-model simulation | a *trace* (`seeded_atlas_trace`) |
