# Phase-atlas pipeline: user guide

The code is in `src/pipelines/seeded_atlas_standalone/`; the one-file entry point is `run_seeded_atlas.m` at the
repository root. Before the first run, `setup_atlas` (section 11) makes the traced `bvp4c` copy from your own MATLAB;
`run_seeded_atlas.m` also does this by itself when the copy is missing or stale.

Words used here (see also the README's "Paper ↔ code" table): a **node** is a grid point of the plane, i.e. one split
of the two conserved quantities between the halves; a **root** is a quasi-steady state found at a node; a **family
label** such as `1+/0` names the sheet (section 5); the **library** is the set of isolated half-cell steady states
used as starting profiles; a **trace** is a full-model simulation started from a root.

The software version (`sa_settings('version')`, also in `VERSION`) is printed in every header and stored in every
run folder. Run folders written by another version still load, plot and trace, but cannot be resumed, and their
traces are recomputed once.

## 1. What it does

The cell is split into two halves, and each conserved quantity (S. pombe: GEF and PAK; M. xanthus: RomR and MglB) is
shared between them. At each node the pipeline looks for roots: solutions of the whole-cell problem of Eq. (8) that
pass every validity check (convergence, positivity, continuity at the middle, the half amounts, residual).

Roots are found in four ways:

- **cold**: a starting profile made of two isolated library halves joined at the middle;
- **continuation**: a step to a neighbouring node from a root already found (the paper's "each solution serves as
  the starting profile at neighbouring points"), halving the step after a failed call;
- **per-root pass**: every distinct root takes one step to each neighbour where no matching root is known (on for
  S. pombe, off for M. xanthus);
- **rescue**: a retry of a node no step reached (default: extrapolation from two roots in a row, "secant").

Every attempt is one `bvp4c` call. A call that returns is not automatically a root: it must pass the checks. Failed
calls are kept, with the reason in plain words.

## 2. Quick start with the wrapper

1. Open `run_seeded_atlas.m` and edit the block at the top: `model`, `par`, `res`, `points_per_axis`, `run_mode`,
   and optionally `library_dir`, `run_dir`, `budget`, `settings_overrides`, `traces`.
2. Run the file. Before any work it prints the model, the parameters, the library folder (exists / will be built),
   the grid (point count or step, the step per axis and the node count), the run folder and what will happen to it,
   the plot path and the trace folders. When the run folder would be refused (`legacy` or `refused`, section 7), it
   stops there, before any library is built.
3. It then builds or reuses the half-cell library, runs (or resumes, or reuses) the run, saves the 3-D plot and runs
   the listed traces. Rerunning the file after a finished run solves nothing; it only adds new traces.

Four things to know:

- **The grid is `points_per_axis`** (shipped: 51): an odd point count per axis, each axis over its own quantity's
  admissible range, so each axis gets its own step. One count per axis is also allowed (`[25 99]`). For a fixed step
  instead, set `points_per_axis = []`; the step is then the model's default (0.02 S. pombe, 0.15 M. xanthus) or
  `settings_overrides.lattice.step`. A run folder holding another grid is refused (`SA:RESUME:E_SETTINGS`), never
  mixed; the default folder name shows the grid, so a run on another grid gets its own folder.
- **Library build.** When the three library files are missing, the wrapper runs Stages 1–3 of
  `run_atlas_surface_pipeline.m` with the script's own calls and arguments (`atlas_build_libraries`). A build opens a
  parallel pool and can take hours. An existing library is reused only when all three files exist and their stored
  parameters equal `par` bit for bit; otherwise the wrapper stops with an error (use another `library_dir`). Write
  `par` with the script's own expressions (`2.38/4`, `2*0.246*7`); the wrapper prints the parameters at full
  precision so a mismatch is visible.
- **`budget` in the edit block is always the limit, in every mode.** The run stops at `max_calls` calls in total or
  after `max_wall_s` seconds of this process, whichever comes first. Running the file again resumes; after a
  `max_calls` stop, raise `max_calls` first. Each field it sets replaces the same field from the mode (smoke 100
  calls / 1 h, pilot 500 calls / 2 h) and from `settings_overrides.budget`; budget fields it does not set, such as
  `settings_overrides.budget.wall_s`, are kept. Production has no budget of its own, so the shipped values
  (500 calls / 7,200 s) apply to a production run too: raise them for a full atlas.
- **The default run folder name shows the mode, the parameters and the grid**:
  `data/standalone_runs/<model>_<mode>_<parameters>_<grid>_<hash>`. Parameters are listed in the order of `par`,
  names without underscores, values printed with `%g` (`.` → `p`, `-` → `m`). The grid part is `n51` for
  `points_per_axis = 51` (`n25x99` for one count per axis) or `step0p02` for a fixed step. For example
  `original_polarity_pilot_a0p595_b0p185_c75p8_d0p552_e0p828_f1_Dm0p05_Dc5_L7_MBG14_MP3p444_n51_92260295` or
  `guzzo_noxa_pilot_L6_rspeed0p35_dbafactor1_rhoR1_rhoB1_dnhalfrange1p8_res25_n51_6d40c23e`. The 8-character hash
  covers the model, the mode and every parameter at full precision, so two parameter sets that print alike still
  get different folders. Nothing that is solved or recorded depends on the name.

**Screen output.**

- **Progress lines.** A `[progress] stage:` line when a stage starts (inputs, resume replay, cold search,
  refinement, continuation, per-root pass, rescue, reports, library build, plot, each trace, compaction) and a
  `[progress] status:` line every `report.progress_interval_s` seconds (default 300; `Inf` = stage lines only), with
  the progress inside the stage where known, calls used of `budget.max_calls`, roots so far and the elapsed time.
- **bvp4c's mesh-limit warning** ("Unable to meet the tolerance without using more than N mesh points ...") is not
  printed. It is recorded for every call (`attempts.csv`, the `convergence_warning` check); a call that hits the limit
  is refused by the convergence checks. `settings_overrides.report.solver_warnings = true` prints it again.

Without the wrapper:

```matlab
setup_atlas;                                        % from the repository folder: puts the package on the path
                                                    % and makes or checks the traced bvp4c copy
s = sa_settings('template', 'guzzo_noxa', 'smoke');
s.run_dir = '/path/to/new/folder';
r = seeded_atlas_run(parameters, struct('directory', library_outputs_folder), s);
```

## 3. Inputs and where defaults come from

- **parameters** (`par`): the model's parameters (docs/MODELS.md). S. pombe: `a b c d e f Dm Dc L M_BG M_P`.
  M. xanthus: `L r_speed dba_factor rho_R rho_B dn_half_range`, plus `res` (the wrapper adds it).
- **library**: the half-cell outputs folder holding `branch_HSS_latest.mat` (no-peak), `branch_HPSS_latest.mat`
  (half-peak) and `branch_MPSS_latest.mat` (full-peak). M. xanthus also accepts a confirmation file
  (`library_confirmation`).
- **settings**: one struct; only `model`, `mode` and `run_dir` are needed.

Each setting is filled in layers, later layers winning: package, model (the model's tested values: grid step
0.02 / 0.15, solver tolerances, starting-mesh size, cold mesh cap 10,000, per-root pass), mode, then yours.
`settings_effective.txt` lists every setting with its value and source (`package`, `model`, `mode`, `user`,
`derived`). Unknown fields and contradictory combinations are refused before any call, and no folder is created.

### Which half shapes are used as starting profiles (`seed_modes`)

Each half of a cold starting profile is one shape from the library. `seed_modes` lists the shapes to try. A shape is
named as the **left** half is drawn, pole to cell middle: `0` (no-peak), `n+` or `n-` for n monotone runs (n = 1, 2,
3, ...; `+` when the half starts high at the pole, `-` when it starts low). The same list is used for the right half,
which is drawn in the opposite direction, so its drawn name is the **mirror** of the shape (the sign flips for odd n
and stays for even n). For example, `1-` on both halves gives `1-/1+`, and `3+` on the right half is drawn `3-`.

- **Default: `{'0','1+','2-'}`** (no-peak, half-peak at the pole, full-peak): the paper's three shapes.
- **`'all'`**: every shape the library folder offers: the three above, `1-` (the half-peak turned around), `2+` (the
  full-peak shifted so its peaks sit at both ends), and both signs of every non-empty higher-mode file (next
  subsection). Each empty mode file that is skipped is named in a printed line. `'all'` cannot be mixed with names.
- Set it with `settings_overrides.seed_modes = {'0','1+','1-','2-','2+'};` or `= 'all';`. A changed list is a new run.
- Families (sheets to seed) are the ordered pairs of the listed shapes, kept when a library cell holds both halves
  and when listed in `families`. The number of pairs grows as the square of the number of shapes (23 shapes: 529).
- **Pre-run table.** Before the first call the run prints one line per seedable family with the expected label, the
  number of library seeds, the cold calls and a rough number of hours, one line per pair without seeds
  (`0 (no cell with both halves)`), then totals. The hours count only cold calls, at the model's measured rate
  (S. pombe about 3.7 s per call at L = 7); read them as an order of magnitude.
- **Refusals (before any call).** An unknown or repeated name, or `'all'` mixed with names (`SA:SEED:E_SHAPE`); a
  missing library file (`SA:SEED:E_FILE`, naming the file and the `max_mode` that builds it); a named mode file whose
  `run_status` is `'empty'` (`SA:SEED:E_EMPTY`); a library file built at other parameters or on another grid
  (`SAP:LIBRARY` or the shared-grid message).
- `1-` and `2+` need no new file; they are made on load from the half-peak and full-peak files.

### Building higher-mode libraries (`n` = 3 and up)

Shapes with more peaks than the paper's three come from extra files `branch_MODE<n>_latest.mat`. A mode-n file is
built from one half-peak piece 1/n of the half long, copied n times, after the three base libraries exist, in the
same folder:

```matlab
atlas_build_libraries('original_polarity', par, res, out_dir, struct('max_mode', 12, 'pool_workers', 1));
```

- **From `run_seeded_atlas.m`:** set `max_mode` in the edit block (default 2 = none). Existing mode files are reused.
  `seed_modes = 'all'` uses only the mode files present; when there are none, the run warns
  (`SA:SEED:W_NO_MODE_FILES`) that only the base shapes are used.
- **`out_dir` must be given** when `max_mode` > 2. `pool_workers` creates the parallel pool with that many workers (a
  pool of another size already open is refused). M. xanthus is built by the same call with its own `par` and `res`.
- **The build stops after the first empty mode** (no cell found a piece); its file is still written with
  `run_status` `'empty'`. "Empty" means the search found nothing, not that the mode does not exist.
- **An unused `MODE<n>` file changes nothing**: it is not loaded unless a listed shape needs it and does not block a
  resume.

### Roots of an earlier run as starting solutions (`import`)

A new run can start from the roots of a finished run folder (the **source**). Each source root goes to the nearest
node of the new grid and is solved there **once**, from its stored solution (no cold cap, no halving, no retry). An
accepted import is a starting root like a cold root: continuation and the per-root pass step from it. The source
must have the same model and parameters; its grid, region, library, solver and search settings may differ (they are
listed in the summary). In the wrapper:

```matlab
import_from = '<a finished run folder>';               % the edit block of run_seeded_atlas.m
settings_overrides.routes.cold = false;                % import only: on the same grid cold would repeat the source
settings_overrides.import.source_edges = 'none';      % or 'same_or_stronger' or 'all' (same grid only; below)
```

- **Mapping.** Nearest node in steps of the new grid (ties to the lowest row); a root farther than
  `import.max_map_steps` (default 1) from every node is skipped (`no_node_near`), one outside `continuation.box` is
  skipped (`outside_box`). The new grid may be the same or finer, nested or not.
- **Mirror images** (`import.mirror`, default true). The left/right mirror of every source root (halves swapped,
  coordinate negated) is also imported at the mirror node, unless its mirror twin is already planned. Both models
  define the reflection; for a model without one, a default `true` becomes `false` with a printed line, and a `true`
  you set is refused (`SA:IMPORT:E_OPTION`).
- **Starting mesh** (`import.guess_mesh`): `'source'` (default: the stored solution on its own mesh) or `'rule'` (the
  continuation rule, `continuation.guess_mesh`).
- **Reusing the source's failed steps** (`import.source_edges`, default `'none'`). Only on the same grid, with the
  same `solver` and `validity` settings, the same `continuation.guess_mesh`, `guess_nodes` and `step_rule` and the
  same model code (else refused, `SA:IMPORT:E_OPTION`; in focus mode the default `'auto'` becomes `'none'` with a
  note). `'same_or_stronger'` skips a continuation step that the source already tried and refused with at least as
  many halvings (the per-root pass skips every source refusal); `'all'` also skips steps that only the source's
  per-root pass refused (weaker evidence). A refused rescue extrapolation is not reused. Skips are counted in the
  summary; a node whose only steps were skipped is `refused` with the reason "refused in the source run (not posed
  here)".
- **The plan** is built without solving, before the run folder is created: every root of the source is read and every
  planned source call file hashed. It is stored in `internal/settings_effective.mat` and never rebuilt on resume. The
  run prints `import: N planned (D direct, M mirror), skipped ...` before the first call.
- **Selecting roots.** Only source roots listed in `import.root_ids` (a cell row such as `{'k:1,4/R000001'}`, or a
  CSV file with a `root_id` column, e.g. rows of the source's `roots_index.csv`) or whose family label is listed in
  `import.families` (`{'1+/0'}`) are imported; both empty (default) imports every root. Mirrors are made for the
  selected roots only. An id or label the source does not hold is refused (`SA:IMPORT:E_SELECTION`).
- **A killed source** (`import.source_stopped = true`). A running source cannot be told from a killed one, so set this
  only when the source will not write again. Its calls without an `end` line in the journal are ignored (with the
  roots they made); a journal written in the last 10 minutes is refused. The source is then identified by the hashes
  of the planned call files, so resuming it later does not block this run. Reuse of its failed steps comes from its
  last checkpoint `queue.mat` only; it needs at least one checkpoint (`attempts.mat`).
- **Refusals (before any folder):** no finished run at `import.source_dir`, an old flat folder, or the run folder
  itself (`SA:IMPORT:E_SOURCE`); a source still running or killed (`E_LIVE`, unless `import.source_stopped` is set; a
  source stopped by its budget is allowed); another model (`E_MODEL`), other parameters (`E_PARAMETERS`), other
  states (`E_STATES`); nothing planned (`E_EMPTY`); the option refusals above (`E_OPTION`); an import in pilot mode.
- **Reports.** `summary.txt` gains an Import section: the source and its differing settings; planned, skipped,
  called, accepted, invalid and unreturned rows; re-solve times; refused imports; label mismatches (information
  only); the largest map distance; root counts per node. A refused import is listed, never retried. The default run
  folder name ends in `_imp<stamp8>`.
- **What an import shows:** that one call converged there from the stored solution, nothing more.

### Focus mode: chosen roots as starting solutions, restricted to their sheets

`run_mode = 'focus'` is import with three differences: you pick the roots, nothing else is used as a start, and the
search is restricted ("fenced") to the sheets you care about. The paper's L = 30 atlas (Fig. 5, S5) was refined this
way.

- **Defaults in focus:** `import_from` is required; `routes.cold = false`; `budget` is required, as in production;
  `import.source_edges = 'auto'` (`'same_or_stronger'` when the grid, solver, validity, step settings and model code
  equal the source's, else `'none'` with a note); `continuation.families = 'selected'` (below). A value you set wins.
  Every other setting is production's.
- **The fence** (`continuation.families`; not the top-level `families`, which selects library pairs for cold
  starts): a list of family labels. A root whose label is outside the list is stored, reported and counted, but no
  continuation, per-root pass or rescue steps from it. A step from inside the fence may land on another sheet; that
  root is kept, only not stepped from. `'selected'` (the focus default) is the labels of the selected source roots,
  plus their mirror labels when `import.mirror` is on; `{}` turns the fence off (with a note). If `'selected'` gives
  no label at all, the run is refused (`SA:IMPORT:E_FENCE`). The fence can also be set in the other modes (default
  `{}`). The summary lists roots stored outside the fence and the steps it withheld.
- **What the fence costs.** A branch that crosses into another sheet and back is not followed through the other
  sheet.
- **No half-cell library.** A focus run neither builds nor reads a library. The wrapper prints
  `library none (mode focus)` and ignores `library_dir`, `res`, `max_mode` and `seed_modes` (M. xanthus also
  `library_confirmation`; its `res` comes from the source's parameters). Two things a run otherwise takes from the
  library come from the source's `internal/settings_effective.mat`: the field scale used to compare profiles (the
  settings it depends on, `lattice.anchor`, `identity.field_scale_floor` and `model_options.hss_branch`, must equal
  the source's, else `SA:IMPORT:E_FOCUS_SOURCE`) and the map from shape names to library families (the source's
  `seed_modes`). No cold targets are computed, `routes.cold = true` is refused (`SA:IMPORT:E_FOCUS_COLD`), and every
  node counts as inside the library range (`search.ring_last`, `rescue.ring` and `budget.ring_calls` have no
  effect). A source that lacks these stored items is refused before any folder (`SA:IMPORT:E_FOCUS_SOURCE`). A
  resume derives both again from the source and refuses on any difference or when the source cannot be read
  (`SA:RESUME:E_FOCUS`); a completed focus run whose source moved still plots and traces, with a note.

Worked example: an earlier production run on 51 points per axis found a sheet `2-/1-` (full-peak/half-peak) worth
refining on a grid twice as fine, without spending calls on the other sheets:

```matlab
run_mode    = 'focus';
points_per_axis = 101;                                  % the finer grid (each root goes to its nearest node)
import_from = '<the earlier run folder>';
budget      = struct('max_calls', 20000, 'max_wall_s', 3*86400);
settings_overrides.import.families = {'2-/1-'};          % or import.root_ids = '<a CSV of chosen roots_index.csv rows>'
% the fence defaults to {'2-/1-', '1+/2-'} (the sheet and its mirror); reuse of failed steps resolves to 'none'
% (another grid), with a printed note
```

## 4. Modes and reading the reports

| mode | cold targets | continuation | budget |
|---|---|---|---|
| smoke | the central nine (the node of equal halves and its 8 neighbours) | inside that 3 × 3 box | 100 calls, 1 h |
| pilot | the central nine plus a 7 × 7 sample over the library's range | along the two grid lines through the centre | 500 calls, 2 h |
| production | one cold call per library seed (`per_seed`, below) | the whole grid | required (the wrapper's `budget`; raise it for a full atlas) |
| focus | none: the selected roots of an earlier run (section 3) | the whole grid, fenced to the selected sheets | required |

**Production cold targets (`per_seed`).** A library seed is a family and one library cell (its mirror cell
follows). Several nodes can map to the same cell; the seed gets one cold call, at the nearest node (ties: the first in
grid order), except that a seed including the centre node is called there. The other nodes are listed as skipped by
policy, not refused. A second root at a skipped node is then looked for only by continuation and the per-root pass;
rescue fills only nodes with no root. With 101 points per axis and a res-50 library this is about a quarter of the
cold calls of `settings_overrides.targets.policy = 'all_seedable'` (one cold call at every seedable pair and node).

**`summary.txt`**: stop reason, calls, targets per family (accepted / invalid / no return), roots, the failing
checks, why calls failed (in plain words, with counts), refinement, identity pairs, the search (node classes, steps,
continuation chains, rescue). Under `per_seed` it prints **MULTIPLICITY NOT TESTED BY COLD**; under any other policy
except `all_seedable`, **MULTIPLICITY NOT TESTED**. Either way a node's root count is a lower bound.

**`pilot_report.txt`** (pilot mode, written at every stop; rebuild with `sa_report('pilot', run_dir)`): 1 settings;
2 refusals by check and route; 3 cost per call and the share of calls that hit the mesh limit (NMax), the number to
watch; 4 mesh along continuation chains; 5 yield per route; 6 the reach map and where positivity or convergence is
lost; 7 identity and multiplicity; 8 the production projection (calls, MATLAB-hours, disk).

> **The projection is an order-of-magnitude guide.** It samples the grid along two lines, so it can miss expensive
> regions or regions without solutions. Set budgets with a margin and watch the run's own summary.

## 5. Results

```
<run_dir>/
  README.txt  summary.txt  roots_index.csv  attempts.csv  settings_effective.txt  [pilot_report.txt]
  plots/      roots_3d.fig, roots_3d.csv, roots_3d_key.mat
              sheets_2d/<family>.fig, all_families.fig, sheets_2d_key.mat
              (relabelled copies: roots_3d_mv-<variable>.*, sheets_2d_mv-<variable>/)
  traces/     <name>/trace.mat, roots_3d_trace.fig, trace_summary.txt, roots_3d_trace_key.mat
  compare/    (after seeded_atlas_compare)
  internal/   the pipeline's state: do not edit, move or delete
```

- **Roots.** `roots_index.csv` has one row per root (`root_id`, node, coordinate, routes, family label).
  `R = seeded_atlas_roots(run_dir)` loads them with every verdict re-derived from the stored calls: use it for any
  analysis, and never load a file in `internal/roots/` directly. It reads the root files one at a time, so memory
  holds one root file plus the list. `struct('sol', true)` also loads each root's solver solution, which holds the
  four sources and so the flux through the middle (Eq. 9).
- **Calls.** `attempts.csv`: one row per call, its checks and, in the last column, the reason in plain words. The
  `file` column is relative to `internal/`.
- **When the files are written.** Each call's own file (`internal/attempts/`) is written as soon as the call
  returns; an accepted call's root file (`internal/roots/`) right after it. During a run, `attempts.csv`,
  `summary.txt` and `roots_index.csv` are rewritten every `report.files_interval_s` seconds (default 600) and at the
  stop. `internal/attempts.mat` and `queue.mat` are checkpoints, written every `checkpoint.every_calls` calls
  (default 25) and at every stop. A resume rebuilds the run from the call files, replaying every stored call, and
  checks each checkpoint and root file against the replay.
- **After a kill.** The three files can be up to one interval old, and a root file can miss the last call's result.
  When the next resume's replay reaches its end, it rewrites such a root file, with two journal lines
  (`note shard ... repair starts`, `note shard ... rewritten`); a `repair starts` line without its partner means the
  process stopped during or just after that write, and the next resume repeats the repair if needed. A kill during a
  write can leave a temporary file `<name>.tmp_<pid>` (or `<name>.mat.tmp_<pid>`) beside the file; the tools never
  read it, and each resume counts these files in its journal and lists them in `internal/resume_<n>.mat`.
- **At the stop** the three files are rendered twice, from the running cache and from scratch; the one from scratch
  is written, with a WARNING line if they differ.
- **Main variable and family labels.** `model_options.main_variable` is the quantity the family label and the plots'
  third axis use: a `+`-separated sum of the model's fields (S. pombe `u q v w p`, default `'u'`; M. xanthus
  `Rm Rc Bm Bc`, default `'Bm'`), e.g. `settings_overrides.model_options.main_variable = 'u+q'`.
- **Family labels (shape names).** A solution is named `<left>/<right>` from the main variable, each half read left
  to right as drawn on the cell (left half: pole → middle; right half: middle → pole):
  - `0` (no-peak): (max − min of the half) / ref is below `model_options.shape_flat_tol`, where ref is the main
    variable's value in the uniform state at the centre node (from the library, interpolated). Defaults: S. pombe
    0.5606, M. xanthus 0.2.
  - `n+` / `n-` otherwise: n monotone runs, where a turn counts only when the profile moves back by at least
    `shape_turn_frac` (0.2) of the half's max − min; `+` when the half's left end is high, `-` when it is low. So a
    peak at the half's left end is `1+`, at its right end `1-`, in its interior `2-`, at both ends `2+`. In the
    paper's words: `1+/0` = half-peak/no-peak, `0/1-` = no-peak/half-peak, `0/0` = no-peak/no-peak, `2-/0` =
    full-peak/no-peak, `2-/1-` = full-peak/half-peak, `1+/1-` = half-peak/half-peak; `1-/1+` is a single peak at the
    cell middle.
  - The cell middle counts as a peak only where the whole cell has a local maximum within `shape_edge_frac` (2 %) of
    a half on either side of it; otherwise a half's last run rising toward the middle is the other half's tail and is
    flattened before the half is named.
  - The mirror image of a root has the halves swapped and the signs of odd modes flipped (`1+/0` ↔ `0/1-`,
    `2-/1-` ↔ `1+/2-`).
  The label is a description: it never chooses a starting profile, a node or an order (every continuation step
  starts from the parent root's own solution). Because the records depend on it, the main variable and the label
  settings are part of a run's identity (changing them is a new run).
- **Plot.** `plots/roots_3d.fig` (open with `openfig`; rotate with the figure toolbar): one colour per family, the two
  plane coordinates plus `dn2` = ∫₀¹ cos(2πη) m dη of the main variable m (S. pombe `d_G`, `d_P`, `dn2`, with the
  view and axis stretch of Fig. 3B; M. xanthus `dn_R`, `dn_B`, `dn2`). The conversion to the paper's axes is in the
  README's "Paper ↔ code" table. `plots/sheets_2d/` holds the same roots in the plane: one `<family>.fig` per family
  and `all_families.fig` with all of them. File names spell `+` as `p`, `-` as `m` and `/` as `_` (`1+/0` →
  `1p_0.fig`). Colours are fixed per label (`1+/0` blue, `0/1-` orange, `0/0` teal, as in Fig. 3A). Coordinates are
  computed from each root's stored profile; `plots/roots_3d.csv` holds them. The plot is a display, not evidence: it
  reads the root files without re-verifying them and makes no solve.
  `seeded_atlas_plot(run_dir)` redraws the plots when the roots or the plot code changed; `struct('verify', true)`
  verifies first.
  **Another main variable, without rerunning:** `seeded_atlas_plot(run_dir, struct('main_variable', 'u+q'))`
  recomputes every root's label and `dn2` with that variable and writes `plots/roots_3d_mv-u+q.*` and
  `plots/sheets_2d_mv-u+q/`; the CSV keeps the run's own label as `stored_family`. Pass the same `main_variable` to
  `seeded_atlas_trace` when a trace's `at` comes from such a figure. Run folders written by older versions show their
  stored `HSS`/`HPSS`/`MPSS` labels; relabel them the same way (`'u'`, M. xanthus `'Bm'`).
- **Traces (full-model simulations).** Started from one root, listed in the wrapper's `traces`:

  | field | meaning | default (S. pombe / M. xanthus) |
  |---|---|---|
  | `name` | folder `traces/<name>/` | required |
  | `root` | a `root_id` from `roots_index.csv` | `root` or `at` required |
  | `at` | a position `[x y z]` in the plot's coordinates (instead of `root`) | |
  | `perturbation` | `struct('kind','none')` or `struct('kind','cos','amplitude',a,'mode',m)`, 0 < a < 1 | none |
  | `t_out` | output times | `[0:0.05:50, 51:1:5000]` s / `linspace(0, 200, 601)` min |
  | `nx` | whole-cell mesh (odd) | 1001 / 241 |
  | `profile_every` | full profiles every Nth frame (first and last always) | 10 |
  | `reltol`, `abstol` | integrator tolerances | 1e-8, 1e-10 / 1e-6, 1e-8 |

  **Starting from a position.** Click a root in `plots/roots_3d.fig`: the data tip shows `at` (ready to paste) and
  `root`. The trace starts from the root closest on the figure (each axis scaled by its plotted range; an axis below
  1 % of the longest counts as 1 %); ties go to the first root in `roots_index.csv`. The chosen root, its distance
  and the next-nearest root are printed and kept in `trace_summary.txt`. If a resume adds a nearer root, the trace is
  integrated again and the reason names both roots.

  The initial condition is the root's own solution on the whole cell; `cos` raises one field of each conserved pair
  and lowers the other by a·min·cos(mπη), so every total is unchanged. S. pombe is integrated with `ode15s` and stops
  when stationary (max|dY/dt| ≤ 1e-8) or when a total drifts by more than 1e-6. M. xanthus is integrated with
  `pdepe` to the last output time; `trace_summary.txt` prints a WARNING when a total drifts by more than 1e-4. Each
  trace folder has `trace.mat` (frames, coordinates, totals, drift, stored profiles), `trace_summary.txt` and
  `roots_3d_trace.fig` (the root plot with the trace: a circle at the start, a square at the end). An unchanged trace
  is reused without integration; a changed entry is integrated again and the old file kept as
  `trace_previous.mat`. Editing comments never invalidates a trace. Sizes: a default S. pombe trace stores about
  20 MB (M. xanthus about 0.5 MB); raise `profile_every` to shrink it.
- **Compaction (optional).** After a completed run, `seeded_atlas_compact(run_dir)` (or `compact_after = true` in the
  wrapper) removes the files of calls that found no root, listed with their SHA-256 in `internal/compaction.txt`. It
  refuses an unfinished run and verifies the roots first. Afterwards the run cannot be resumed (`E_COMPACTED`), and
  the tools that count every call (`sa_report('cutset', ...)`, the pilot report rebuild) refuse; roots, plots, traces,
  comparisons and import keep working. A run that **stopped cleanly without completing** (on disk space, budget or
  wall time) is compacted only with `seeded_atlas_compact(run_dir, struct('stopped_ok', true))`; a killed run is
  refused even then (resume it to a stop first). On large runs most of the disk is in the files of calls that found
  no root.

## 6. Tuning

A changed scientific setting is a **new run** in a new folder. Set only what you change, in the wrapper's
`settings_overrides` or on the settings struct.

| Symptom | Setting | Notes |
|---|---|---|
| Grid: a point count per axis or a fixed step | `lattice.points_per_axis` (odd n ≥ 3, scalar or one per axis; the wrapper sets it from `points_per_axis`) | each axis gets its own step from its admissible range: the open range (−R, R) with the outermost node half a step inside, or [−cap, cap] when `lattice.max_abs_coordinate` is set. The centre node must exist, so the count is odd (use 49 or 51, not 50). Not the library's `res`, which spans 95 % of the range and may be even. Not together with `lattice.step` |
| Grid too coarse or too expensive | `points_per_axis`; or `lattice.step` with `points_per_axis = []` | default step 0.02 (S. pombe), 0.15 (M. xanthus) |
| Only part of the plane is wanted | `lattice.max_abs_coordinate` (scalar or one per axis) | keeps nodes with \|coordinate\| ≤ cap; must keep the centre |
| Many calls hit the mesh limit (high NMax share, `convergence_*` refusals) | `solver.nmax` and `solver.mesh_reject` (≤ nmax); `continuation.guess_nodes` | S. pombe: 802 starting-mesh nodes is the first lever (default 401); M. xanthus default 482 |
| Cold solves that need a large mesh are cut | `cold.mesh_cap` (cold calls only; default 10,000; `Inf` = off; ≤ `solver.mesh_reject`) | a call cut by the cap is not evidence that no root exists; `sa_report('cutset', run_dir, c)` counts, without solving, which cold calls a cap c cuts |
| Coexisting roots missed, or the per-root pass too expensive | `root_pass.enabled`; `root_pass.max_calls_per_round`, `root_pass.max_rounds` (default 10000: the pass runs until a round adds nothing) | M. xanthus: off by default; `max_calls_per_round` ends the pass for the current search pass (a cost cap) |
| Continuation chains cut | `continuation.max_lineages`, `continuation.max_calls_per_lineage`, `multiseed.max_roots_per_node` (all default 10000, i.e. effectively unbounded) | lower them to bound the cost |
| Start from an earlier run's roots | wrapper `import_from`; `import.mirror`, `import.guess_mesh`, `import.max_map_steps`, `import.source_edges` | section 3; `routes.cold = false` for import only |
| Only some roots of the earlier run are wanted | `import.root_ids`, `import.families`; `run_mode = 'focus'` | section 3 |
| Calls spent on sheets you do not care about | `continuation.families` (the fence; focus default `'selected'`; `{}` = none) | roots outside the fence are kept, never stepped from |
| The earlier run was killed | `import.source_stopped = true` | only when it will not write again |
| Nodes next to reached ones stay empty | `rescue.rungs = {'diagonal', 'secant', 'deep'}` (that order); `rescue.ring = true` (rescue also outside the library's range) | default secant only, no rescue outside the library's range |
| Rescue spends calls and rescues nothing | `routes.rescue = false`, or `rescue.max_calls_per_node` | |
| Projection above your hours or disk | `lattice.step`, `lattice.points_per_axis`, `lattice.max_abs_coordinate`; `targets.policy = 'central_nine'` | `central_nine` maps only the part connected to the centre and does not test multiplicity |
| Too many cold targets of one family | `targets.thin = struct('families', {{'HSS_HSS'}}, 'stride', 2)` with `targets.policy = 'all_seedable'` | multiplicity then only partly tested; refused with `per_seed` |
| Cold target choice | `targets.policy`: `per_seed` (production default), `all_seedable`, `central_nine`, `pilot_sample`, `explicit` (with `targets.list`), `nearest_anchor` (with `targets.per_family`) | only `all_seedable` gives a full cold multiplicity census |
| The run must stop cleanly in pieces, or run longer | `budget.max_calls`, `budget.max_wall_s` and `budget.wall_s`; the wrapper's `budget` sets the first two | stops at whichever comes first; may change on resume |
| Roots too loose or too tight | `solver.reltol`, `solver.abstol` | abstol below reltol; the refinement tolerances must stay tighter |
| Progress lines too frequent or too rare | `report.progress_interval_s` (`Inf` = stage lines only) | default 300 s; may change on resume |
| Report files refreshed too often or too rarely | `report.files_interval_s` (`Inf` = only at the stop) | default 600 s; may change on resume |
| Checkpoints written too often or too rarely | `checkpoint.every_calls` | default 25; may change on resume |
| Debugging a solve | `report.solver_warnings = true` | prints bvp4c's warnings |
| Which half shapes start the search | `seed_modes` (default `{'0','1+','2-'}`; `'all'`; section 3) | read the pre-run table first |
| Cold calls ordered from the centre outward | `targets.policy = 'per_seed_centre'` | `per_seed`'s targets, nearest the centre first |
| Long cells | `scaling.reference_length` (S. pombe 14 µm, M. xanthus 6 µm) and the counts below | see "The floor rule"; a value you set always wins |
| Labels or `dn2` should use another variable | `model_options.main_variable` (`'u'`, `'q'`, `'u+q'`, ...; M. xanthus `'Bm'`, `'Rm+Bm'`, ...) | a new run; for an existing run relabel the plot instead (section 5) |
| Too many or too few halves named no-peak, or small wiggles counted as turns | `model_options.shape_flat_tol` (S. pombe 0.5606, M. xanthus 0.2), `model_options.shape_turn_frac` (0.2) | a new run |

There is no setting for a search region away from the centre: moving `lattice.anchor` is a diagnostic only.

**The floor rule (point counts grow with L).** The default mesh and node counts were set at L ≤ 14 µm (S. pombe) and
L ≤ 6 µm (M. xanthus), the *reference length*. A longer cell holds longer profiles, so above it these counts are
scaled by L / reference length: `model_options.guess_nodes`, `continuation.guess_nodes`, `model_options.shape_nodes`
(nodes: 1 + ceil((N − 1) L / L_ref)) and `solver.nmax`, `solver.mesh_reject`, `cold.mesh_cap` (ceil(N L / L_ref)).
Their source in `settings_effective.txt` becomes `derived`. At or below the reference length nothing changes. A
value you set yourself is never overwritten. S. pombe at L = 30 gets 859 starting-mesh nodes, `nmax` and
`mesh_reject` 42,858 and a cold cap of 21,429. Not scaled: `identity.n_comparison_nodes`, `validity.n_dense`,
`identity.relsup_nodes`. A mode-n shape uses 1/n of the nodes per piece; when (`guess_nodes` − 1) / n is below 50 the
run warns (`SA:SEED:W_PIECE_NODES`). A call file at `nmax` 42,858 is about 19 MB.

## 7. Stops, reuse and resume

Stop reasons: `completed`, `budget_calls`, `budget_wall`, `budget_ring`, `disk_floor` (free space below
`checkpoint.min_free_gib`, default 8 GiB). A killed process can also be resumed.

The wrapper decides from the run folder (`sa_resume('status', ...)`, which writes nothing). A `legacy` or `refused`
folder stops the wrapper before the library step:

| status | when | what the wrapper does |
|---|---|---|
| `new` | the folder is missing or empty | runs |
| `legacy` | a folder in an old flat format (no `internal/`) | refuses; the folder still loads, plots and traces into another folder (`opts.out`) |
| `refused` | the parameters, a scientific setting, the software, MATLAB, model, library or code version differ, the run stopped on an integrity error, or the folder is not empty and holds no run | refuses, naming every difference, with zero calls and the folder untouched: choose a new `run_dir`. When the parameters differ, the library loader's message (`SAP:LIBRARY` for S. pombe, `SAG:LIBRARY` for M. xanthus) may appear instead |
| `completed` | the run finished (also after compaction) | no solve and no replay; plot and traces only |
| `budget_exhausted` | stopped at `budget.max_calls` and the budget is not raised | no solve: raise `budget.max_calls` to continue |
| `resumable` | any other stop | resumes: every stored call is replayed and checked without solving, then the run continues |

Only `budget.*`, `checkpoint.*`, `report.*` and `import.source_dir` (a moved source) may differ on resume. Everything
else, including the mode and the grid, is scientific. Edits to the plot, trace, wrapper, compaction, loader or
`post/` code never block a resume.

A finished run is checked like a resume. When package code the run depends on (the `sa_*.m` files,
`seeded_atlas_run.m`, the model module, `solver/` and the private folders) changed since the run started, even a
`completed` folder is `refused` (`E_CODE`), and the wrapper stops before the plot and the traces. To plot and trace
such a folder, call the tools directly; they make no solve and write only `plots/` and `traces/<name>/`:

```matlab
setup_atlas;                                   % from the repository folder
seeded_atlas_plot(run_dir);                    % plots/roots_3d.fig, .csv, _key.mat and plots/sheets_2d/
seeded_atlas_trace(run_dir, struct('name', 'C1', 'root', 'k:1,4/R000001'));   % traces/C1/
```

Through the wrapper instead, add `settings_overrides.resume = true` and
`settings_overrides.resume_code_override = '<why the code changed>'` with the run's own settings: the status is then
`completed`, and still nothing is solved. A call cut off by a kill leaves an `interrupted` record that counts toward
`budget.max_calls`; the resume runs it again.

**An import run and its source.** The source is identified by a stamp (SHA-256 over its `settings_effective.mat`,
`attempts.mat`, `queue_final.mat`, `journal.txt`, `compaction.mat` and the call-file count). A resume of the new run
recomputes it and refuses a changed source (`E_IMPORT`); each import call also checks the SHA-256 of its source call
file. A moved source resumes with `import.source_dir` set to its new place. Resuming or compacting the source after
the new run started changes the stamp and so blocks the new run's resume; a `completed` new run stays `completed`
(with a note), so it can still be plotted and traced. Plotting or tracing the source does not change the stamp. A run
from a killed source (`import.source_stopped`) is identified by the hashes of its planned call files instead.

**A focus run** stores in `internal/settings_effective.mat` (`focus_inputs`) where its field scale and shape-name map
came from (the source path and the SHA-256 of its `settings_effective.mat`, the import identity, the field scale, the
source's `seed_modes`). A resume derives them again and refuses on any difference or when the source cannot be read
(`E_FOCUS`); a `completed` focus run stays `completed`, with a note.

**The per-root pass and `root_pass_cap`.** `root_pass.max_calls_per_round` is a cost cap: when reached, the pass ends
for the current search pass (`root_pass_cap`); a later pass (when rescue added a node) may pose the rest.

Directly: `s.resume = true; s.budget.max_calls = 1000; seeded_atlas_run(p, manifest, s)`. For production after a
pilot, `[s, dropped] = sa_settings('from_run', pilot_dir, 'production')` keeps your own pilot settings and drops the
mode's (listed in `dropped`); set `s.run_dir` and the budgets.

## 8. Comparing with a reference

```matlab
C = seeded_atlas_compare(run_dir, reference);   % writes run_dir/compare/<name>.csv and <name>_map.txt
```

`reference` is a struct (or a `.mat` file holding one) with one entry per node: coordinate, profile on its own mesh,
usable flag and label (see `help seeded_atlas_compare`). Each node is classed agree (≤ 1e-6), band, different
(≥ 1e-2), ours-only or reference-only. A reference never starts or steers a run; the classes are a report, not a
pass/fail test. The converters `'from_myxo_sheet'` and `'from_stageb_registry'` read the authors' own reference files,
which are not included.

## 9. What a result does not mean

A root is a quasi-steady state that passed the checks at these tolerances. It is not shown to be stable, two
neighbouring roots are not shown to lie on one branch, and two roots count as one when their profiles differ by less
than 1e-3 (as in the paper). A node without a root is not proof that none exists, and a call cut by the cold mesh cap
is not evidence either. A trace's end state is not shown to be an attractor. The third plot coordinate, `dn2`, is a
display coordinate; it is not shown to separate every pair of sheets.

## 10. Tested settings and timing

The settings at which the pipeline was tested, and those used in the paper, are listed in [MODELS](MODELS.md)
("Tested settings"). Another L, other parameters, a new grid or another MATLAB release is untested: start with a
pilot run. Timing notes for long cells (from runs at L = 30, library grid `res = 15`):

- **One pool worker gives no parallelism but doubles the MATLAB hours** (the client and the worker both run): the
  base library build at L = 30 took about 1 hour per process, about 2 MATLAB-hours. Budget hours as processes times
  wall time.
- **The half-peak library covers little of the plane at L = 30, `res = 15`** (16 of 225 cells; no-peak 225, full-peak
  70). Sheets with a half-peak half are then thinly started, and "no root" for them says little.
- **The scaled cold mesh cap cuts most higher-mode cold calls at L = 30** (74 of 81 at the cap 21,429). If those
  roots are wanted, raise `cold.mesh_cap` (`Inf` = no cap; it must not exceed `solver.mesh_reject`).

## 11. Operating notes

- **Generated solver copy.** The pipeline calls a traced copy of your MATLAB's own `bvp4c` that records every
  proposed mesh (`solver/sa_bvp4c_traced.m` and `solver/private/`, with the record `solver/sa_bvp4c_generated.txt`).
  These files are not in git and no MathWorks file is shipped. Make them with `setup_atlas` (from the repository
  folder); `run_seeded_atlas.m` calls `setup_atlas('quiet')` before every run and remakes the copy when it is missing
  or does not match your MATLAB. `setup_atlas('check')` verifies only and prints one line, e.g.
  `setup_atlas: OK  R2025b (25.2.0.3150157) maca64  bvp4c 597839cef613  toolbox=reference  deval=replica
  selfcheck=pass(7/7)`, or `setup_atlas: FAIL <check>: <reason>`.
- **Other MATLAB releases.** A release is accepted when its `bvp4c.m` can be patched and the copy passes the
  self-check: on small fixed problems it gives exactly the same results as your untouched `bvp4c`, its mesh hooks fire
  where expected, its text equals your `bvp4c.m` except the renamed function line and two added lines, and the
  warning and error identifiers the pipeline uses exist. R2025b Update 4 is labelled `reference`, any other release
  `self-checked`. A `.p` or MEX `bvp4c`, or a rewritten one, is refused with a named check (`E_TOOLBOX`, `E_ANCHOR`,
  `E_CLOSURE`, `E_SELFCHECK`). The fast `deval` replica is used only with the reference `deval.m`
  (`deval=replica`); otherwise your own `deval` is used (`deval=toolbox`). The self-check tests the patch, not the
  numerics: another release may find other roots, call counts or meshes.
- **MATLAB updates and other machines.** Every run records the MATLAB version (with the Update), the architecture
  and the toolbox hashes. A run cannot be resumed after a MATLAB Update or on another architecture (`E_VERSION`) or
  another `bvp4c` ("the MATLAB toolbox changed since this run started"). Finished folders still load, plot, trace
  and compare.
- **Do not run the tests while a run is in progress in the same repository folder.** The tests that make traced
  copies (`sa_t61p_traced_solver.m`, `sa_t102_release_selfcheck.m`) work only in scratch folders under
  `data/experiments/`; the note is a precaution.
- The half-cell script `run_atlas_surface_pipeline.m` always reruns its own simulation stage (`force_trace = true`);
  the wrapper never runs that stage.
