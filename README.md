# polarity-atlas

MATLAB pipelines that map the steady states of two cell-polarity models (fission yeast, `original_polarity`, and the
myxobacterial Guzzo_noXA model, `guzzo_noxa`) on a cell split into two halves:

- the **uncoupled atlas** (`run_atlas_surface_pipeline.m`) builds libraries of one-half-cell steady states (flat,
  half-peak, middle-peak and higher modes) over a grid of half-cell pool masses;
- the **seeded coupled atlas** (`run_seeded_atlas.m`) uses those libraries as starting guesses ("seeds") to find
  steady states of the whole two-half cell at every node of a lattice of pool splits, with continuation, rescue,
  resume after a stop or a kill, plots and full-PDE traces.

Version `1.0.0-rc1` (see `VERSION`, `CHANGELOG.md`). Shared privately for testing; see `LICENSE`.

## Requirements

- MATLAB with the Java VM (never start MATLAB with `-nojvm`).
- **Tested on R2025b Update 4, macOS arm64. Other releases are self-checked, not tested; the minimum release is
  unverified** (by reading function introduction dates: R2020a for everything, R2019a for library builds, coupled
  runs, plots and traces).
- Toolboxes: Optimization (`fsolve`), Signal Processing (`findpeaks`) and Parallel Computing (`parfor`,
  `parfeval`) are needed **only** by the fission-yeast library build and the 12-stage uncoupled script. The Guzzo_noXA
  path, coupled runs, plots and traces need base MATLAB only.
- **Wen-Han and JB: please run `setup_atlas('check')` (after `setup_atlas`) and send the one line it prints**, for
  example
  `setup_atlas: OK  R2025b (25.2.0.3150157) maca64  bvp4c 597839cef613  toolbox=reference  deval=replica  selfcheck=pass(7/7)`
  or `setup_atlas: FAIL <check>: <reason>`. This line is the first evidence on another MATLAB.

## Setup

From the package folder in MATLAB:

```matlab
setup_atlas            % puts src/pipelines/seeded_atlas_standalone on the path and makes the traced bvp4c copy
setup_atlas('check')   % verifies only; prints the one-line verdict
```

`setup_atlas` makes `src/pipelines/seeded_atlas_standalone/solver/sa_bvp4c_traced.m` (and `solver/private/`) from
**your own MATLAB's** `bvp4c.m`: the same file with a renamed function line and two added lines that record the mesh
sizes it proposes. It then runs a self-check (about 12 tiny `bvp4c` calls, a few seconds). These files are git-ignored;
no MathWorks file is in this repository. `run_seeded_atlas.m` calls `setup_atlas('quiet')` itself before every run.

## Quick start

1. **A small coupled run.** Open `run_seeded_atlas.m`, edit the block at the top and run the file. For a first try:
   `run_mode = 'smoke';` and a small library grid, e.g. fission yeast `res = 15;` (or the Guzzo_noXA example in the
   block, with `res = 11;`). When the library folder is empty the script builds the library itself (Stages 1-3 of the
   uncoupled pipeline; the fission-yeast build opens a parallel pool and can take an hour or more). It prints every
   path before any work.
2. **Plot.** The run saves `plots/roots_3d.fig` and the 2-D sheets; redraw with `seeded_atlas_plot(run_dir)`.
3. **One trace.** A full-PDE integration from one root:
   `seeded_atlas_trace(run_dir, struct('name', 'T1', 'root', '<a root_id from roots_index.csv>'))`, or list it in the
   script's `traces` block.
4. **A library built directly.** `setup_atlas` adds only the coupled package, so first add the uncoupled pipeline:
   ```matlab
   addpath('src/solvers/atlas_surface_pipeline_dependencies');
   % Guzzo_noXA also: addpath('src/solvers/atlas_surface_pipeline_dependencies/models/guzzo_noxa'); addpath('src/models/guzzo_noxa');
   m = atlas_build_libraries('guzzo_noxa', par, 25, fullfile(pwd, 'data', 'my_library'));
   ```
5. **The 12-stage uncoupled pipeline:** edit and run `run_atlas_surface_pipeline.m` (see `docs/UNCOUPLED_GUIDE.md`).

Documents: `docs/USER_GUIDE.md` (the coupled pipeline; settings, modes, results, resume), `docs/UNCOUPLED_GUIDE.md`
(the uncoupled pipeline), `docs/MODELS.md` (models, parameters, commissioned scope), `docs/SPEC.md` (the product
contract and claim ceiling).

## Folder map

```
setup_atlas.m                    path + traced bvp4c copy + self-check
run_seeded_atlas.m               the coupled pipeline (edit the block at the top)
run_atlas_surface_pipeline.m     the 12-stage uncoupled pipeline
src/pipelines/seeded_atlas_standalone/              the coupled package (models/, post/, solver/, private/)
src/solvers/atlas_surface_pipeline_dependencies/    the uncoupled pipeline (incl. models/guzzo_noxa/)
src/models/guzzo_noxa/                               the shared Guzzo_noXA core
tests/seeded_atlas_standalone/   the portable tests: run_portable_tests
docs/                            USER_GUIDE, UNCOUPLED_GUIDE, MODELS, SPEC
data/                            created by the first run; git-ignored; never shipped
```

**Outputs** go to `data/` in the package folder (libraries `data/<parameters>/outputs`, coupled runs
`data/standalone_runs/<name>`), unless you give a folder. `data/` is git-ignored.

**Tests:** `addpath('tests/seeded_atlas_standalone'); run_portable_tests` (about 15 minutes; stub models and
correctors, plus about 140 tiny `bvp4c` calls in self-checks of scratch copies). It prints PASS, FAIL or SKIP per
test; scratch goes to `data/experiments/portable_tests_<stamp>/`.

## Release policy (other MATLAB releases)

The pipeline accepts any MATLAB release whose `bvp4c.m` can be patched and whose patched copy passes the self-check;
it refuses only with a message naming the failed check.

**Guaranteed on a release that passes the self-check:**
- the traced solver is that release's own `bvp4c` plus two observation lines, shown by the text diff and by
  bit-identical outputs on the self-check problems;
- the mesh hooks fire where the pipeline's per-call consistency rule expects;
- the warning and error identifiers the pipeline relies on exist;
- the pipeline code, settings, acceptance gates and records are the same;
- `deval` is the package's replica only on the reference `deval.m`; otherwise it is your own (parity by construction);
- no run mixes releases, Updates or architectures;
- every run names the release, the architecture and the toolbox hashes.

**Not guaranteed:**
- the same roots, call counts, meshes, timings or library coverage as R2025b;
- bit-identity across releases, or across machines on one release;
- the commissioned scope, which holds at R2025b Update 4, maca64 only;
- that the package works at all on another release: **untested**. The self-check line is the evidence a collaborator
  produces.

The uncoupled pipeline uses your own `bvp4c`, `pdepe`, `fsolve` with no self-check. Its libraries are new inputs on
each release, kept apart by their SHA-256 in the coupled run's identity.

## Platform notes

1. **Line endings.** `.gitattributes` (`* -text`) keeps every file byte-exact on every platform; code hashes and the
   manifest depend on it. Do not re-save files with an editor that converts line endings.
2. **Windows (untested).** Drive paths (`C:\...`) count as relative where the wrapper changes folder, so the run stays
   in your current folder; this is harmless and recorded in the run folder. Default run-folder names are about 90
   characters, and paths under `data/standalone_runs/<name>/roots/` can exceed 260 characters: clone to a short path or
   enable long paths. Test T82 skips its macOS-only case.
3. **Never start MATLAB with `-nojvm`** (hashing and file handling use Java).
4. **`metadata_log.csv`** appears in the current folder after a fission-yeast library build (it is git-ignored).
5. **Toolboxes:** see Requirements.
6. **Pool size is your choice:** open `parpool(n)` first or pass `opts.pool_workers`; the default pool uses all
   cores. One worker gives no parallelism but costs two MATLAB processes.
7. **Running tests during a run:** do not run the shipped tests while a run is in flight in the same package folder.
   The self-check tests (T61p, T102) generate only into scratch and never move or rewrite the installed copy.
8. **Release and platform:** tested on R2025b Update 4, macOS arm64; other releases are self-checked, not tested; the
   minimum release is unverified. Please report the `setup_atlas('check')` line.
9. **MATLAB Updates:** a run cannot be resumed after a MATLAB Update. Finish it first, or start a new `run_dir`.

## Known limitations

- **Uncoupled HPSS libraries can hold unconverged entries** (development review R17, A-1): the HPSS sweeper accepts
  the solver's return without a convergence gate, so some stored HPSS cells have `sol.stats.maxres` above the solver's
  RelTol (in the development libraries: L = 7, 12 of 1,096; L = 7.94, 30 of 1,142; L = 14, 198 of 866). As coupled
  seeds this is harmless (a failed cold call); existence read from the uncoupled masks is overstated.
- **The MPSS / mode-n constructor's slope gates test only flipped pieces** (R17, A-2): `seam_slope` and
  `trough_slope` are read at the BVP ends, where zero flux makes them exactly 0, so the interior join's true slope is
  never checked. Legacy behaviour, kept unchanged.
- Windows is untested; other MATLAB releases are self-checked, not tested.
- No resume across a MATLAB Update or another architecture.
- An unpatchable `bvp4c` (a `.p` or MEX file, or a rewrite that moved its anchor lines) is refused; there is no
  mode without the mesh trace.

## Provenance and notes

Exported from the development repository (private), commit `5d05db64` of `main`, on 2026-10-31; history not
included. `EXPORT_MANIFEST.txt` lists every file's SHA-256.

Comments in the code cite documents of the development repository (packets, reports, `docs/...`, `figures/...`);
those documents are not included. The user guide keeps its section numbers so that the code's references to
"USER_GUIDE section N" still point to the right place (`docs/USER_GUIDE.md`).

`src/pipelines/seeded_atlas_standalone/post/guzzo_noxa/Guzzo_noXA.m` is a byte copy of Wen-Han's script, kept because
it enters the trace key; it is never executed.

**Credits: to be completed before public release.**
