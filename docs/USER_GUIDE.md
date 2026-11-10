# Standalone seeded coupled atlas: user guide

The code is `src/pipelines/seeded_atlas_standalone/` (software 1.0.0-rc2). The one-file entry point is
`run_seeded_atlas.m` in the package folder. Before the first run, `setup_atlas` (section 11) makes the traced
`bvp4c` copy from your own MATLAB; `run_seeded_atlas.m` also does this by itself when the copy is missing or stale.

The software version is the text returned by `sa_settings('version')` (in `sa_settings.m`); it is printed in every
header and stored in every run folder. It is not a setting and users do not change it. A new version means: run
folders written by an older version still load, plot and trace, but cannot be resumed, and their traces are
recomputed once (the version is part of every trace key).

## 1. What it does

A cell is split into two halves, and each conserved pool (fission yeast: P and G; Guzzo_noXA: R and B) is shared
between them. A **lattice node** is one way of splitting the pools; its coordinate is the left-minus-right difference.
At each node the pipeline looks for **roots**: steady states of the coupled two-half problem that pass every validity
check (convergence, positivity, seam match, pool masses, residual).

Roots are found from the uncoupled seed library in four ways:

- **cold**: a guess assembled from two library halves;
- **continuation**: a step to a neighbouring node from a root already found, halving the step after a failed call;
- **per-root pass**: every distinct root takes one step to each neighbour where no matching root is known (on for
  fission yeast, off for Guzzo_noXA);
- **rescue**: a retry of a node no walk reached (default: secant extrapolation from two roots in a row).

Every attempt is one `bvp4c` call. A call that returns is not automatically a root: it must pass the checks. Failed
calls are kept, with the reason in plain words.

## 2. Quick start with the wrapper

1. Open `run_seeded_atlas.m` and edit the block at the top: `model`, `par`, `res`, `points_per_axis`, `run_mode`,
   and optionally `library_dir`, `run_dir`, `budget`, `settings_overrides`, `traces`.
2. Run the file in MATLAB. It prints the model, the parameters, the library folder (exists / will be built), the
   coupled grid (point count or step, the derived per-axis steps and the node count), the run folder and what will
   happen to it, the plot path and the trace folders, before any work. When the run folder
   would be refused (`legacy` or `refused`, section 7), it stops right there, before any library is built.
3. It then builds or reuses the uncoupled library, runs (or resumes, or reuses) the coupled run, saves the 3-D root
   plot and runs the listed traces. Rerunning the file after a finished run solves nothing; it only adds new traces.

Four things to know:

- **The coupled grid is `points_per_axis`** (shipped: 51), like the uncoupled `res`: an odd point count per pool
  axis, each axis over its own pool's admissible range, so each axis gets its own step. One count per pool is also
  allowed (`[25 99]`). For a fixed step instead, set `points_per_axis = []`; the step is then the model's default
  (0.02 fission yeast, 0.15 Guzzo_noXA) or `settings_overrides.lattice.step`. The grid is a scientific setting:
  a run folder holding another grid is refused (`SA:RESUME:E_SETTINGS`, naming `lattice.points_per_axis`), never
  mixed. The default folder name shows the grid (next point), so a run on another grid gets its own folder.
- **The library build path is argument-verified only.** When the three library files are missing, the wrapper runs
  Stages 1–3 of `run_atlas_surface_pipeline.m` with exactly the script's calls and arguments (`atlas_build_libraries`);
  this equivalence was checked argument for argument, not by rebuilding a library. A build opens a parallel pool
  (`p.use_parallel = true`, as in the script) and can take hours. An existing library is reused only when all three
  files exist and their stored parameters equal `par` bit for bit; otherwise the wrapper stops with an error (use
  another `library_dir`). Write `par` with the script's own expressions (`2.38/4`, `2*0.246*7`); the wrapper prints
  the parameters at full precision so a mismatch is visible.
- **`budget` in the edit block is always the limit, in every mode, production included.** The run stops at
  `max_calls` calls in total or after `max_wall_s` seconds of this process, whichever comes first. Running the file
  again resumes; after a `max_calls` stop, raise `max_calls` first (otherwise nothing is solved). Each field it sets
  replaces the same field from the mode (smoke 100 calls / 1 h, pilot 500 calls / 2 h) and from
  `settings_overrides.budget`; budget fields it does not set, such as `settings_overrides.budget.wall_s`, are kept.
  Production has no budget of its own, so the shipped values (500 calls / 7,200 s) apply to a production run too:
  raise them for a full atlas.
- **The default run folder name shows the mode, the parameters and the grid**:
  `data/standalone_runs/<model>_<mode>_<parameters>_<grid>_<hash>`. The parameter part lists every parameter in the
  order of `par`, under its name without underscores, with the value printed as the uncoupled folder names print it
  (`%g`, `.` → `p`, `-` → `m`). For fission yeast this is exactly the uncoupled folder's parameter text. The grid
  part is the grid the run uses, printed the same way: `n51` for `points_per_axis = 51` (`n25x99` for one count per
  pool), or, with `points_per_axis = []`, the fixed step: `step0p02` (`step0p01x0p05` for one step per pool). It is
  taken after `settings_overrides` (and `from_pilot`), so a grid set there shows too. For example
  `original_polarity_pilot_a0p595_b0p185_c75p8_d0p552_e0p828_f1_Dm0p05_Dc5_L7_MBG14_MP3p444_n51_92260295` (the shipped
  block), `original_polarity_pilot_a0p595_b0p185_c75p8_d0p552_e0p828_f1_Dm0p05_Dc5_L7_MBG14_MP3p444_step0p02_92260295`
  (the fixed step); for Guzzo_noXA,
  `guzzo_noxa_pilot_L6_rspeed0p35_dbafactor1_rhoR1_rhoB1_dnhalfrange1p8_res25_n51_6d40c23e` (the uncoupled Guzzo_noXA
  folders use shorter tags). The 8-character hash covers the model, the mode and every parameter at full precision
  (not the grid), so two parameter sets that print alike still get different folders, and the pilot and the
  production run of one parameter set never share a folder. A folder from an older version has no grid part in its
  name, so the default name never points at it. Nothing that is solved or recorded depends on the name.
  `data/standalone_runs/` is not tracked by git.

**Screen output.**

- **Progress lines.** A `[progress] stage:` line when a stage starts (inputs, resume replay, cold search, refinement,
  continuation, per-root pass and rescue for each phase and search pass, final reports, pilot report, library build,
  plot, each trace, compaction) and a `[progress] status:` line every `report.progress_interval_s` seconds (default
  300; `Inf` = stage lines only). Each line gives the stage, the progress inside it where it is known (cold targets
  called, box nodes holding a root, edges), calls used of `budget.max_calls`, roots so far and the elapsed time.
  Status lines are checked after every call and, for fission-yeast traces, during the integration. They change
  nothing that is solved.
- **bvp4c's mesh-limit warning** ("Unable to meet the tolerance without using more than N mesh points ...") is not
  printed. It stays enabled and is recorded for every call (`attempts.csv`, the `convergence_warning` check); a call
  that hits the limit is refused by the convergence checks as before. A one-line note says so the first time.
  `settings_overrides.report.solver_warnings = true` prints it again; anything else a solve prints is always shown.
  Both settings may change on resume.

Without the wrapper:

```matlab
setup_atlas;                                        % from the package folder: puts the package on the path and
                                                    % makes or checks the traced bvp4c copy; nothing else is added
s = sa_settings('template', 'guzzo_noxa', 'smoke');
s.run_dir = '/path/to/new/folder';
r = seeded_atlas_run(parameters, struct('directory', library_outputs_folder), s);
```

## 3. Inputs and where defaults come from

- **parameters** (`par`): the model's dimensional parameters. Fission yeast: `a b c d e f Dm Dc L M_BG M_P`.
  Guzzo_noXA: `L r_speed dba_factor rho_R rho_B dn_half_range`, plus `res` (the wrapper adds it).
- **library**: the uncoupled outputs folder holding `branch_HSS_latest.mat`, `branch_HPSS_latest.mat` and
  `branch_MPSS_latest.mat`. Guzzo_noXA also accepts a pre-seeding confirmation file (`library_confirmation`).
- **settings**: one struct; only `model`, `mode` and `run_dir` are needed.

Each setting is filled in layers, later layers winning: package, model (that model's commissioned values: lattice
step 0.02 / 0.15, solver tolerances, guess mesh, cold mesh cap 10,000, per-root pass), mode, then yours.
`settings_effective.txt` lists every setting with its value and source (`package`, `model`, `mode`, `user`,
`derived`). Unknown fields and contradictory combinations are refused before any call, and no folder is created.

### Which half shapes are seeded (`seed_modes`, since 0.12)

Each half of a cold guess is one shape taken from the library. `seed_modes` is the list of shapes you want tried.
A shape is named as the **left** half is drawn, pole to cell middle: `0` (flat), `n+` or `n-` for n peaks (n = 1, 2,
3, ...; `+` when the half starts high at the pole, `-` when it starts low). The same list is used for the right half;
the right half is drawn in the opposite direction, so its drawn name is the **mirror** of the shape (the sign flips
for odd n and stays for even n). For example, seeding `1-` on both halves gives the guess `1-/1+`, and `3+` on the
right half is drawn `3-`.

- **Default: `{'0','1+','2-'}`** (flat, half-peak at the pole, middle peak). This is exactly the seeding of 0.11; an
  unchanged run gives bit-identical results.
- **`'all'`**: every shape the library folder offers: the three above, the new `1-` (the half-peak turned around) and
  `2+` (the middle peak rotated, so its peaks sit at both ends), and both signs of every mode-n file that is not
  empty (next subsection). Each empty mode file that is skipped is named in a printed line. `'all'` cannot be mixed
  with names.
- Set it in the wrapper with `settings_overrides.seed_modes = {'0','1+','1-','2-','2+'};` or `= 'all';`. It is a
  scientific setting (a new run).
- Families are the ordered pairs of the listed shapes, filtered as before by the census (a pair needs a library
  cell with both halves) and by `families`. There is no pairing restriction, so the number of pairs grows as the
  square of the number of shapes (23 shapes: 529 pairs).
- **Pre-run table.** Before the first call the run prints one line per seedable family with the expected drawn
  label, the number of library seeds, the cold calls this run will make and a rough number of hours, and one line per
  pair with no seeds (`0 (no cell with both halves)`), then totals. The hours are only the cold calls at the model's
  measured cold rate (fission yeast 3.7 s per call at L = 7); continuation, refinement and the per-root pass are not
  counted, and a family with a new shape says `no measured rate; legacy rate used`. Read it as an order of magnitude;
  if it is too large, shorten `seed_modes` or set `families` before the first solve.
- **Refusals (before any call).** An unknown name, a repeated name, or `'all'` mixed with names
  (`SA:SEED:E_SHAPE`); a shape whose library file is missing (`SA:SEED:E_FILE`; the message names the file and the
  `max_mode` that builds it); a mode file you name explicitly whose `run_status` is `'empty'` (`SA:SEED:E_EMPTY`); a
  library file built at other parameters or on another grid (`SAP:LIBRARY` or the shared-grid message).
- `1-` and `2+` need no new file; they are made on load from the HPSS and MPSS files.

### Building mode libraries (`n` = 3 and up)

The three base libraries hold flat, half-peak and middle-peak halves. Shapes with more peaks come from extra files
`branch_MODE<n>_latest.mat`. A mode-n file is built by the uncoupled builder from one half-peak piece that is 1/n of
the half long, copied n times; it is built after the three base libraries exist, in the same folder:

```matlab
atlas_build_libraries('original_polarity', par, res, out_dir, struct('max_mode', 12, 'pool_workers', 1));
```

- **From `run_seeded_atlas.m`:** set `max_mode` in the edit block (default 2 = none; for example 12). The wrapper
  passes it to this call for the run's library folder; existing mode files are reused. `seed_modes = 'all'` uses only
  the mode files present: when the folder holds none with n ≥ 3, the run prints a warning
  (`SA:SEED:W_NO_MODE_FILES`) naming the folder and saying that only the base shapes are seeded.
- **`out_dir` must be given** (an empty one is refused when `max_mode` > 2). `max_mode` is the highest n to try
  (default 2 = build nothing new); a mode whose file already exists is reused. `pool_workers` creates the parallel
  pool with that many workers (a pool of another size already open is refused). Guzzo_noXA is built by the same call
  with its own `par` and `res`.
- **The build stops after the first empty mode** (no cell found a piece); its file is still written with
  `run_status` `'empty'`. This is a cost rule: it saves the calls for higher modes, and "empty" is a search result
  (the build did not find a piece), not proof that the mode does not exist.
- **An unused `MODE<n>` file changes nothing**: it is not loaded unless an active shape needs it, it does not enter
  the library fingerprint, and a run folder with such a file present still resumes.

### Roots of an earlier run as seeds (`import`, since 1.0.0-rc2)

A new run can start from the roots of a finished run folder (the **source**). Each source root goes to the nearest
node of the new grid and is solved there **once**, from its stored solution (route, kind and phase `import`; the
run's solver settings, no cold cap, no halving, no retry). An accepted import is a starting root like a cold root:
it may start or join a continuation lineage (under the default `continuation.starts = 'anchor_and_cold'`) and the
per-root pass steps from it in round 1. The source must have the same model and parameters; its grid, box, library,
solver and search settings may differ (they are listed in the summary). Set it in the wrapper:

```matlab
import_from = '<a finished run folder>';               % the edit block of run_seeded_atlas.m
settings_overrides.routes.cold = false;                % import only: on the same grid cold would repeat the source
settings_overrides.import.source_edges = 'none';      % or 'same_or_stronger' or 'all': your per-run choice, same grid only (below)
```

- **Mapping.** Nearest node by the distance in new-lattice steps (ties to the lowest lattice row); a root farther
  than `import.max_map_steps` (default 1) from every node is skipped (`no_node_near`), one outside
  `continuation.box` / region is skipped (`outside_box`). The grid may be the same or finer and not nested.
- **Mirror images** (`import.mirror`, default true). The left/right mirror image of every source root is also
  imported at the mirror node (the halves swapped, the coordinate negated), unless the source's own mirror twin is
  already a planned import (a twin within the identity tolerance at the node nearest −c). Both shipped models declare
  the reflection (module member `reflection`); for a model without one a default `true` becomes `false` with a printed
  line, and a `true` you set is refused (`SA:IMPORT:E_OPTION`).
- **Guess mesh** (`import.guess_mesh`): `'source'` (default: the stored solution on its own mesh) or `'rule'` (the
  continuation guess-mesh rule, `continuation.guess_mesh`).
- **Source-edge knowledge** (`import.source_edges`, default `'none'`; a per-run choice). On the same grid with the
  same `solver` and `validity` groups, the same `continuation.guess_mesh`, `guess_nodes` and `step_rule` and the same
  model version (the corrector and solver code: a refusal depends on it) (else refused, `SA:IMPORT:E_OPTION`; in
  focus mode the default `'auto'` becomes `'none'` with a printed note instead), the run skips edges the source already tried and refused, from the new root its
  imported parent landed on: `'same_or_stronger'` skips a lineage edge only when the source refused it with at least
  as many halvings (h ≥ `continuation.max_halvings`; the per-root pass skips every source refusal); `'all'` also
  lets a lineage skip edges only the source's per-root pass (h = 1) refused (weaker evidence: deeper halving might
  have completed some). A completed source edge is used by the pass only when its child was imported too. A
  refused secant-rescue call carries no knowledge: its guess is extrapolated from a lineage, not taken from the
  parent root's own solution, so it says nothing about a step from that root. Skips are
  counted by category in the summary; a node whose only edges were skipped is `refused` with the reason "refused in
  the source run (not posed here)".
- **The plan** is built with zero solves before the run folder is created (a refusal or a kill during it leaves no
  folder): reading every root of the source (`seeded_atlas_roots`, unverified) and hashing every planned source
  attempt file. It is stored in `internal/settings_effective.mat` (`import`) and never rebuilt on resume. The run
  prints `import: N planned (D direct, M mirror), skipped ...` before the first call.
- **Selecting roots** (`import.root_ids`, `import.families`). Only the source roots whose id is listed in
  `import.root_ids` (a cell row such as `{'k:1,4/R000001'}`, or the path of a CSV file with a `root_id` column, for
  example rows cut from the source's `roots_index.csv`) or whose observed family label is listed in
  `import.families` (`{'1+/0'}`) are imported; both empty (the default) imports every root. Mirror rows are made for
  the selected roots only (a mirror whose twin is not selected is called). An id or label the source does not hold
  is refused (`SA:IMPORT:E_SELECTION`, listing them). The summary gives the lists and the selected / skipped counts.
- **A killed source** (`import.source_stopped = true`). There is no lock or process file, so a running source cannot
  be told from a killed one: set it only when the source will not write again. Its calls without an `end` line in
  the journal are ignored (and the roots they made); a journal written in the last 10 minutes is refused. Its
  identity is the SHA-256 of each planned row's attempt file (attempt files are never rewritten), not the folder
  stamp, so the source may be resumed later without blocking this run's resume (and the default run folder name uses
  the source's `settings_effective.mat` instead, so it stays the same too). Source-edge knowledge comes from its last
  checkpoint `queue.mat` only. The summary says that the source was killed, at which call, and which checkpoint was
  used. The source needs at least one checkpoint (`attempts.mat`).
- **Refusals (before any folder):** no finished run at `import.source_dir`, a flat (0.5 or earlier) folder or the
  run folder itself (`SA:IMPORT:E_SOURCE`); a source still running or killed: a journal start without an end, or a
  `queue_final.mat` count that differs from the attempt files (`E_LIVE`, unless `import.source_stopped` is set; the
  message says so; a budget-stopped source is allowed); another
  model (`E_MODEL`), other parameters (`E_PARAMETERS`, naming the fields), other states (`E_STATES`); import only
  with nothing planned (`E_EMPTY`); the option refusals above (`E_OPTION`); an import in pilot mode (the pilot
  projection has no import model).
- **Reports.** `summary.txt` gains an Import section: source path, stamp and versions, the differing settings,
  planned / skipped / called / accepted (new root or matched) / invalid / no-return for direct and mirror rows, the
  re-solve wall time by returned-mesh bin, the refused imports by root and label, label mismatches (information
  only: labels are descriptions), the largest map distance, and the root counts per node of the source and of this
  run. An import that is refused is listed, never retried, never a stop. The default run folder name ends in
  `_imp<stamp8>`. With no cold roots, refinement is not run.
- **What an import shows:** that one call converged there from the stored solution, nothing more; the imports are
  not shown to be all the roots near those nodes.

### Focus mode: chosen roots as seeds, with a family fence

`run_mode = 'focus'` is import with three differences: you pick the roots, nothing else is seeded, and the search is
fenced to the families you care about.

- **Defaults in focus:** `import_from` (`import.source_dir`) is required (refused without); `routes.cold = false`
  (the imports are the only seeds); `budget` is required, as in production; `import.source_edges = 'auto'`
  (`'same_or_stronger'` when the grid, solver, validity, step-rule settings and model version equal the source's,
  else `'none'` with a printed note); `continuation.families = 'selected'` (the fence, below). A value you set wins.
  Every other setting is production's.
- **The family fence** (`continuation.families`; not the top-level `families`, which names library families for
  cold seeding). A list of observed family labels. A root whose label (the label of the call that created it) is
  outside the list is stored, reported and counted, but it never becomes a lineage's expandable root, never starts a
  lineage (the anchor lineage starts from the first root inside the fence), never starts the per-root pass and is
  never a rescue parent or secant basis. A step aimed from a fenced family may land on another family; that root is
  kept (it is a result), only not stepped from. `'selected'` (the focus default) is the labels of the selected source
  roots plus their mirror labels when `import.mirror` is on; `{}` turns the fence off (focus then prints a note that it
  runs without a fence). If `'selected'` would give no label at all (selected roots without an observed label), the
  run is refused (`SA:IMPORT:E_FENCE`) rather than run unfenced. The fence can also be set in the other modes (default
  `{}` there). The summary gets a fence section: roots stored outside the fence (by label), the lineage expansions,
  lineage starts and per-root pass starts it withheld, and the rescue nodes whose ladder had no parent or basis
  inside the fence (recorded `no_basis`, counted apart from genuine `no_basis` nodes).
- **What the fence costs.** A branch that crosses into another family and back is not followed through the other
  family; roots beyond it are found only if some fenced root reaches them directly.

Worked example: an earlier production run on 51 points per axis found a family `2-/1-` worth refining, and you want
it on a grid twice as fine, without spending calls on the other families:

```matlab
run_mode    = 'focus';
points_per_axis = 101;                                  % the new, finer grid (import maps each root to its nearest node)
import_from = '<the earlier run folder>';
budget      = struct('max_calls', 20000, 'max_wall_s', 3*86400);
settings_overrides.import.families = {'2-/1-'};          % or import.root_ids = '<a CSV of chosen roots_index.csv rows>'
% the fence defaults to {'2-/1-', '1+/2-'} (the family and its mirror); source edges resolve to 'none'
% (another grid), with a printed note
```

The run plans the selected roots (zero solves), prints the selection and the fence, re-solves each root once at its
nearest node of the 101-point grid, and continues only from roots labelled `2-/1-` or `1+/2-`.

## 4. Modes and reading the reports

| mode | cold targets | continuation | budget |
|---|---|---|---|
| smoke | the central nine (the node nearest equal halves and its 8 neighbours) | inside that 3 × 3 box | 100 calls, 1 h |
| pilot | the central nine plus a 7 × 7 sample over the library's range | along the two lattice lines through the centre | 500 calls, 2 h |
| production | one cold call per library seed (`per_seed`, below) | the whole lattice | none of its own: required (the wrapper's edit-block `budget` supplies it; raise it for a full atlas) |
| focus | none: the selected roots of an earlier run are the seeds (`import_from` required; section 3) | the whole lattice, fenced to the selected families | required, as production |

**Production cold targets (`per_seed`).** A library seed is a family and one library cell (its mirror cell follows).
Several grid nodes can map to the same cell; the seed gets one cold call, at the seedable node nearest the cell
(ties: the first in lattice order). The centre node (the anchor) is the exception: a seed that includes it is called
there, so the centre gets the same cold calls as before. The other nodes of the seed are listed as policy-skipped,
not refused. Only the cold calls change. A second root at one of a seed's skipped nodes is then looked for only by
continuation (capped lineages) and the per-root pass (neighbours of roots found); rescue only fills nodes with no root.
Continuation lineages other than the anchor's also start from different cold roots. With
101 points per axis and a res-50 library this is about a quarter of the cold calls. For one cold call at every
seedable (family, node) pair, as before, set `settings_overrides.targets.policy = 'all_seedable'`.

**`summary.txt`**: stop reason, calls, targets per family (accepted / invalid / no return), roots, the failing
checks, **why calls failed** in plain words with counts, refinement, identity pairs, the search (node classes,
edges, lineages, rescue). Under `per_seed` it gives the seed count and the seedable pairs not called cold, and
prints **MULTIPLICITY NOT TESTED BY COLD**: other roots at a seed's other nodes are looked for only by the later
phases. For any other policy except `all_seedable` it prints **MULTIPLICITY NOT TESTED**: the roots found are those
connected to the cold-seeded nodes. Either way a node's root count is a lower bound.

**`pilot_report.txt`** (pilot mode, written at every stop; rebuild with `sa_report('pilot', run_dir)`): 1 settings;
2 refusals by check and route; 3 cost per call, and the share of calls that hit the mesh limit (NMax), the number to
watch; 4 mesh along continuation chains; 5 yield per route; 6 the reach map and the frontier (where positivity or
convergence is lost); 7 identity and multiplicity; 8 the production projection (calls, MATLAB-hours, disk).

> **The projection may be inaccurate; treat it as an order-of-magnitude guide.** It samples the lattice along two
> lines, so it can miss expensive regions or regions without valid states. Set budgets with a margin and watch the
> run's own summary.

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

- **Roots.** `roots_index.csv` has one row per root (`root_id`, node, coordinate, routes, observed family).
  `R = seeded_atlas_roots(run_dir)` loads them the verified way (every verdict re-derived from the stored calls):
  use it for any analysis. Never load a file in `internal/roots/` directly. It reads the root files one at a time
  and keeps only each root's row and stored profile, so memory holds one root file plus the list, not the whole
  run (the root files of a production run do not fit in memory together). `struct('sol', true)` also loads every
  root's solver solution into the list, which does grow with the run; `seeded_atlas_compare` reads them one at a
  time instead.
- **Calls.** `attempts.csv`: one row per call, its checks and, in the last column, the reason in plain words. The
  `file` column is relative to `internal/`.
- **When the files are written.** Every call's own file (`internal/attempts/`) is written as soon as the call
  returns; when the call is accepted, its node's root file (`internal/roots/`) is written right after it. While a
  run is going, `attempts.csv`, `summary.txt` and `roots_index.csv` are rewritten every `report.files_interval_s`
  seconds (default 600, checked after each live call) and always when the run stops. `internal/attempts.mat` and
  `queue.mat` are checkpoints, written every `checkpoint.every_calls` calls (default 25) and at every stop, but
  never during a resume's replay (a stop before the replay has finished keeps the stored ones). A resume does not
  start from them: it rebuilds the run from the call files, replaying every stored call, and compares each
  checkpoint and each root file with the replay at the call count it was written at.
- **After a kill.** The three files can be up to one interval old, and a root file can lack the discovery of the
  call the kill followed (a new root, or one more discovery of a root it already holds). When the next resume's
  replay reaches its end (not when it stops earlier, for example on its budget), it rewrites such a root file, with
  two journal lines per file: `note shard ... repair starts` before the write and `note shard ... rewritten` after
  it. A `repair starts` line without its `rewritten` line means the process stopped during that write, or just
  after it and before the `rewritten` line; the next resume repeats the repair if the file is still behind, and
  writes nothing if the write had finished. The three files are rewritten at the first live call after
  the replay when the replay took longer than `report.files_interval_s`, otherwise at the first interval after it,
  or at the stop. A kill during a write can leave a temporary file `<name>.tmp_<pid>` beside the file it was
  replacing (for a `.mat` file too: `<name>.mat.tmp_<pid>`, as MATLAB R2025b saves it). The resume itself never
  reads it, nor do `seeded_atlas_roots`, `seeded_atlas_plot`, `seeded_atlas_compact` and the pilot report; each
  resume counts those in `internal/` and its `attempts/` and `roots/` folders in its journal line
  `note resume <n>: ... stray temporary file(s)` and lists them in `internal/resume_<n>.mat`.
- **At the stop** the pipeline renders the three files twice, from its running cache and from scratch, and writes
  the one from scratch with a WARNING line (screen and journal) if they differ.
- **Main variable and family labels.** `model_options.main_variable` is the quantity the family label and the
  plots' third axis look at: a `+`-separated sum of the model's fields, spaces allowed (fission yeast `u q v w p`,
  default `'u'`; Guzzo_noXA `Rm Rc Bm Bc`, default `'Bm'`), e.g. `settings_overrides.model_options.main_variable =
  'u+q'`. Anything else is refused with the allowed field names.
- **Mode names (the family label, since 0.11).** A call's solution is named `<left>/<right>` from the main
  variable, each half read left to right as drawn on the full cell (left half: pole → cell middle; right half:
  cell middle → pole):
  - `0`: flat, (max − min of the half) / ref below `model_options.shape_flat_tol`, where ref is the main
    variable's value in the uniform (HSS) state at the anchor: the library HSS interpolated bilinearly to the
    anchor (stored in the run's field scale, `anchor_hss_values`). Defaults: fission yeast 0.5606 (= 0.25 /
    0.44595, the earlier absolute 0.25 in u at the L = 7 setting); Guzzo_noXA 0.2 (depleted Bm halves whose
    whole amplitude is a few per cent of the uniform Bm read as flat; the earlier relative contrast called a
    half flat almost never).
  - `n+` / `n−` otherwise: n monotone runs, where a turn counts only when the profile moves back by at least
    `shape_turn_frac` (0.2) of the half's max − min; `+` when the half's left end is high, `−` when it is low
    (written `-`). So a peak at the half's left end is `1+`, at its right end `1-`, in its middle `2-`, at both
    ends `2+`, 1.5 periods `3+` / `3-`. Examples: `1+/0` (the earlier HPSS_HSS), `0/1-` (HSS_HPSS), `0/0`
    (HSS_HSS), `2-/0` (MPSS_HSS), `1+/2-` (HPSS_MPSS), a single peak at the cell middle `1-/1+`.
  - The cell middle counts as a peak only where the whole cell has a local maximum within `shape_edge_frac`
    (2 %) of a half on either side of it (a top at most a quarter of a turn higher just beyond that window
    still counts, so a broad middle peak a few per cent off the middle reads `1-/1+`); otherwise a half's last run rising toward the middle is the
    neighbour's tail, not a peak, and is left out (flattened to its lowest value) before the half is named.
  - The mirror image of a root has the halves swapped and odd modes' signs flipped (`1+/0` ↔ `0/1-`,
    `2-/1-` ↔ `1+/2-`).
  The label sets each call's `observed_family` and, in the per-root continuation pass, the unit's family
  string (a name maps to a library family: `0` → HSS, a pole peak of mode 1 → HPSS, `2-` → MPSS; any other name
  uses the family of the call that found the root), which enters the unit keys a resume replays and the new
  root's recorded seed family. It never chooses a guess, a node, an order or a cap: every continuation step
  starts from the parent root's own solution. Because records and resume identity depend on it, the main
  variable and the label settings are scientific settings (a new run).
- **Plot.** `plots/roots_3d.fig` (open with `openfig` or a double-click; rotate with the figure toolbar): one colour
  per observed family, the two lattice coordinates plus Δn₂ = ∫₀¹ cos(2πη) m dη of the main variable m (fission yeast
  ΔG, ΔP, Δn₂ with Fig 3b's view and stretch; Guzzo_noXA Δn_R, Δn_B, Δn₂). The axis label names the variable, e.g.
  Δn₂ (u). Markers are filled and half transparent, so overlapping sheets show. `plots/sheets_2d/` holds the same
  roots in the mass-imbalance plane (the 3-D plot's two lattice axes): one `<family>.fig` per observed family and
  `all_families.fig` with every family overlaid and a legend, same colours. Sheet file names spell a mode name
  with `p` for `+`, `m` for `-` and `_` for `/` (`1+/0` → `1p_0.fig`, `2-/1-` → `2m_1m.fig`). Colours are fixed per
  name; a mode name keeps the colour of the earlier family it replaces (`1+/0` the HPSS_HSS blue, `0/1-` the
  HSS_HPSS orange, `0/0` the HSS_HSS teal, ...); any other name gets a colour fixed by the name itself. Coordinates are computed from each
  root's stored profile; `plots/roots_3d.csv` holds them and the gap to the node's lattice coordinate. The plot
  reads the root files decoded but not re-verified: it is a display, not evidence. It makes no solve and reads the
  root files one at a time, once for the 3-D plot and the sheets together, so a production run plots in about the
  memory of one root file plus the plotted roots.
  `seeded_atlas_plot(run_dir)` redraws the 3-D plot and the sheets when the roots or the plot code changed
  (`plots/roots_3d_key.mat`, `plots/sheets_2d/sheets_2d_key.mat`; a missing file is redrawn too);
  `struct('verify', true)` verifies first.
  **Another main variable, without rerunning:** `seeded_atlas_plot(run_dir, struct('main_variable', 'u+q'))`
  recomputes every root's family (the current rule, on its stored profile) and Δn₂ with that variable and writes
  `plots/roots_3d_mv-u+q.*` and `plots/sheets_2d_mv-u+q/`, titled "relabelled for plotting; run labels as stored";
  the CSV keeps the run's label as `stored_family`. The run's own files and labels are unchanged. The relabel reads
  the stored profile (on the comparison mesh), not each call's solver solution as the run does; on the folders
  checked the two agreed on every root. Pass the same `main_variable` in `seeded_atlas_trace`'s third argument when
  a trace's `at` comes from a relabelled figure.
  **Folders written before 0.11** show their stored HSS / HPSS / MPSS labels (a line printed by the plot says
  so). To see them with mode names, relabel: `struct('main_variable', 'u')` (Guzzo_noXA `'Bm'`); the uniform-state
  reference at the anchor is recomputed from the folder's library files (nothing is solved); if they cannot be
  read, the stored HSS of the native cell nearest the anchor is used and a line says so. Folders
  written before 0.10 also record no main variable: their stored labels were computed on u / Bm and Δn₂ is
  computed on u / Bm (fission yeast figures before 0.10 used u + q; the printed line says so).
- **Traces.** A full-PDE integration started from one root, listed in the wrapper's `traces`:

  | field | meaning | default (fission yeast / Guzzo_noXA) |
  |---|---|---|
  | `name` | folder `traces/<name>/` | required |
  | `root` | a `root_id` from `roots_index.csv` | `root` or `at` required |
  | `at` | a position `[x y z]` in the plot's coordinates (instead of `root`) | |
  | `perturbation` | `struct('kind','none')` or `struct('kind','cos','amplitude',a,'mode',m)`, 0 < a < 1 | none |
  | `t_out` | output times | `[0:0.05:50, 51:1:5000]` / `linspace(0, 200, 601)` min |
  | `nx` | full-domain mesh (odd) | 1001 / 241 |
  | `profile_every` | full profiles every Nth frame (first and last always) | 10 |
  | `reltol`, `abstol` | integrator tolerances | 1e-8, 1e-10 / 1e-6, 1e-8 |

  **Starting from a position.** Click a root in `plots/roots_3d.fig`: the data tip shows `at` (its coordinates as
  `[x y z]`, ready to paste into the entry's `at`) and `root` (its `root_id`). The trace starts from the root closest
  on the figure: each axis difference is divided by that axis's plotted range (largest minus smallest coordinate
  over the run's roots; an axis that is nearly flat on the figure, below 1 % of the longest axis, counts as that
  1 %), then the usual distance is taken; ties go to the first root in `roots_index.csv`. The
  chosen root, its coordinates and the distance (scaled, and raw per axis), with the next-nearest root, are printed
  and kept in `trace_summary.txt`. The trace is reused while the entry and the chosen root stay the same; when a
  resume adds roots and another root becomes the nearest, the trace is integrated again and the reason names the old
  and the new root.

  The initial condition is the root's own solver solution assembled on the full domain; `cos` raises one field of
  each mass-conserving pair and lowers the other by a·min·cos(mπη), so every pool is unchanged. Fission yeast
  integrates with `ode15s` and stops when stationary (max|dY/dt| ≤ 1e-8) or when a pool mass drifts by more than
  1e-6. Guzzo_noXA integrates with `pdepe` to the horizon; `trace_summary.txt` prints a WARNING when a pool mass
  drifts by more than 1e-4. Each trace folder has `trace.mat` (frames, coordinates, masses, drift, stored profiles),
  `trace_summary.txt` and `roots_3d_trace.fig` (the root plot with the trace: a circle at the start, a square at the
  end). A trace whose entry, root, settings and code are unchanged is reused without integration; a changed entry is
  integrated again and the old file kept as `trace_previous.mat`. Editing comments in the code never invalidates a
  trace; `roots_3d_trace_key.mat` decides when the figure is redrawn (after a resume added roots, for example).
  Sizes: a default fission-yeast trace stores about 20 MB (Guzzo_noXA about 0.5 MB); raise `profile_every` to
  shrink it. A root plot is tens of kB; a trace figure adds about 25 kB per 1,000 trace points.
- **Compaction (optional).** By default every call's file is kept, so a run can be resumed and every root
  re-verified. After a completed run `seeded_atlas_compact(run_dir)` (or `compact_after = true` in the wrapper)
  removes the files of calls that found no root, listed with their SHA-256 in `internal/compaction.txt`. It refuses an
  unfinished run and verifies the roots before removing anything. Afterwards the run cannot be resumed, and the
  tools that count every call (the cut-set recount, the pilot report rebuild) refuse; roots, plots, traces and
  comparisons keep working.
  A run that **stopped cleanly without completing** (a `disk_floor`, `budget_calls` or wall-time stop: its
  `queue_final.mat` is present and counts every attempt file) is compacted only on request,
  `seeded_atlas_compact(run_dir, struct('stopped_ok', true))`; without it the refusal says so. It prints that the run
  can no longer be resumed, and `compaction.mat` / `compaction.txt` record that it was compacted from a stopped run,
  with the stop reason. A killed run (no final queue, a count that differs, or a call started without an end) is
  refused even then: resume it to a stop first. A compacted stopped run still plots, traces and serves as an import
  source (section 3; a focus run can start from it), but cannot be resumed (`E_COMPACTED`). On large runs most of the
  disk is in the files of calls that found no root, so this frees most of the space.

## 6. Tuning

A changed scientific setting is a **new run** in a new folder. Set only what you change, in the wrapper's
`settings_overrides` or on the settings struct.

| Symptom | Setting | Notes |
|---|---|---|
| Grid: a point count per axis (the wrapper's `points_per_axis`, shipped 51) or a fixed step | `lattice.points_per_axis` (an odd n ≥ 3, scalar or one per pool; the wrapper sets it from `points_per_axis`) | Each axis gets its own step from that pool's admissible range: the open range (−R, R) with the outermost node half a step inside, or the closed range [−cap, cap] when `lattice.max_abs_coordinate` is set. The anchor is the middle node, so the count must be odd (50 is refused: use 49 or 51). This is not the uncoupled `res`, which is a closed grid over 95 % of the range and allows even counts. Not together with `lattice.step`, and not with a moved anchor |
| Lattice too coarse or too expensive | `points_per_axis` in the wrapper; or a fixed `lattice.step` (scalar or one per pool) with `points_per_axis = []` | the package default without the wrapper's count: 0.02 (fission yeast), 0.15 (Guzzo_noXA) |
| Only part of the lattice is wanted | `lattice.max_abs_coordinate` (scalar or one per pool) | keeps nodes with \|coordinate\| ≤ cap; the cap must keep the anchor |
| Many calls hit the mesh limit (NMax share high, `convergence_*` refusals) | `solver.nmax` and `solver.mesh_reject` (≤ nmax); `continuation.guess_nodes` | fission yeast: 802 guess nodes is the documented lever (default 401); Guzzo_noXA default 482 |
| Cold roots that need a large mesh are cut | `cold.mesh_cap` (cold calls only; default 10,000; `Inf` = off; ≤ `solver.mesh_reject`) | a call cut by the cap is a refusal of the capped search, not evidence that no root exists; `sa_report('cutset', run_dir, c)` counts, without solving, which cold calls a cap c cuts |
| Coexisting roots missed, or the pass too expensive | `root_pass.enabled`; `root_pass.max_calls_per_round`, `root_pass.max_rounds` (default 10000 since 1.0.0-rc2, was 8: the pass runs until a round adds nothing) | Guzzo_noXA: opt in with `true`; `max_calls_per_round` ends the pass for the current search pass (a user cost cap) |
| Lineages blocked or cut | `continuation.max_lineages`, `continuation.max_calls_per_lineage`, `multiseed.max_roots_per_node` (all default 10000 since 1.0.0-rc2; were 16, 100 and 4) | the old values cut the L = 30 run's lineages (16 of 16 slots, 10 at the call cap, 575 of 1,286 nodes closed to lineages); set the old values to reproduce a 1.0.0-rc1 run; a pilot's projection then prints the cap column as "unbounded" |
| Start from an earlier run's roots | wrapper `import_from`; `import.mirror`, `import.guess_mesh`, `import.max_map_steps`, `import.source_edges` | section 3; `routes.cold = false` for import only |
| Only some roots of the earlier run are wanted | `import.root_ids` (a cell row of root ids or a CSV path with a `root_id` column), `import.families` (labels); `run_mode = 'focus'` | section 3; an unknown id or label is refused, listing them |
| Calls spent on families you do not care about | `continuation.families` (the fence: a list of labels; focus default `'selected'`; `{}` = none) | roots outside the fence are kept and counted, never stepped from |
| The earlier run was killed | `import.source_stopped = true` | only when it will not write again; its calls without an end line are ignored; a journal written in the last 10 minutes is refused |
| Nodes next to reached ones stay empty | `rescue.rungs = {'diagonal', 'secant', 'deep'}` (that order); `rescue.ring = true` (rescue also outside the library's range) | default secant only, no ring rescue |
| Rescue spends calls and rescues nothing | `routes.rescue = false`, or `rescue.max_calls_per_node` | |
| Projection above your hours or disk | `lattice.step`, `lattice.points_per_axis`, `lattice.max_abs_coordinate`; `targets.policy = 'central_nine'` | `central_nine` maps only the branch connected to the centre and does not test multiplicity |
| Too many cold targets of one family | `targets.thin` = `struct('families', {{'HSS_HSS'}}, 'stride', 2)` with `targets.policy = 'all_seedable'` | keeps a (family, node) pair on a stride; multiplicity is then only partially tested; refused with `per_seed` |
| Cold target choice | `targets.policy`: `per_seed` (production default: one cold call per library seed), `all_seedable` (every seedable pair), `central_nine`, `pilot_sample`, `explicit` (with `targets.list`), `nearest_anchor` (with `targets.per_family`; nearest in lattice steps) | only `all_seedable` gives a cold multiplicity census; `per_seed` leaves it to the later phases |
| The run must stop cleanly in pieces, or run longer | `budget.max_calls` (calls in total), `budget.max_wall_s` and `budget.wall_s` (seconds of this process); the wrapper's edit-block `budget` sets the first two and always applies, in every mode | the run stops at whichever comes first; budgets may change on resume; after a `max_calls` stop, raise it to continue; production has no budget of its own, so raise these for a full atlas |
| Roots too loose or too tight | `solver.reltol`, `solver.abstol` | abstol below reltol; the refinement tolerances must stay tighter |
| Progress lines too frequent or too rare | `report.progress_interval_s` (seconds; `Inf` = stage lines only) | default 300; may change on resume |
| `attempts.csv`, `summary.txt`, `roots_index.csv` refreshed too often or too rarely during a run | `report.files_interval_s` (seconds; `Inf` = only when the run stops) | default 600; they are always written at the stop; changes nothing solved or recorded; may change on resume |
| The checkpoints (`attempts.mat`, `queue.mat`) written too often (slow on a long run) or too rarely | `checkpoint.every_calls` | default 25; a resume replays every stored call either way; may change on resume |
| Debugging a solve: see bvp4c's warnings | `report.solver_warnings = true` | default false: the mesh-limit warning is recorded per call, not printed; may change on resume |
| Which half shapes are seeded; more or fewer cold families | `seed_modes` (default `{'0','1+','2-'}`; `'all'`; section 3) | a new run; read the pre-run table first (families grow as the square of the number of shapes) |
| Cold calls ordered from the centre outward | `targets.policy = 'per_seed_centre'` | per_seed's targets, the cold calls sorted nearest d = 0 first (ties keep the per_seed order); a run option, never a default; the report treats it as `per_seed` |
| Large cells (above the reference length) | `scaling.reference_length` (fission yeast 14 µm, Guzzo_noXA 6 µm; model layer) and the leaves below | see "The floor rule" below; a value you set always wins |
| Families or Δn₂ should look at another variable | `model_options.main_variable` (`'u'`, `'q'`, `'u+q'`, ...; Guzzo_noXA `'Bm'`, `'Rm+Bm'`, ...) | a new run; for an existing run relabel the plot instead (section 5); the flatness reference follows the variable (its uniform value) |
| Too many or too few halves named flat (`0`), or small wiggles counted as turns | `model_options.shape_flat_tol` (fission yeast 0.5606, Guzzo_noXA 0.2; relative to the main variable's uniform value), `model_options.shape_turn_frac` (0.2) | a new run; `shape_amp_tol` and `shape_contrast` (before 0.11) are no longer settings and are refused |

There is no setting for a box away from the centre: moving `lattice.anchor` is a diagnostic only.

**The floor rule (point counts grow with L, since 0.12).** The commissioned point counts were set at L ≤ 14 (fission
yeast) and L ≤ 6 (Guzzo_noXA), the *reference length*. Above it a cell holds longer peaks, so these leaves are
scaled up by L / reference length: `model_options.guess_nodes`, `continuation.guess_nodes`,
`model_options.shape_nodes` (nodes: 1 + ceil((N - 1) L / L_ref)) and `solver.nmax`, `solver.mesh_reject`,
`cold.mesh_cap` (mesh limits: ceil(N L / L_ref)). Their source in `settings_effective.txt` becomes `derived`.
At or below the reference length nothing changes (a run at L = 7 or 14 is bit-identical to 0.11). A value you set
yourself, including `cold.mesh_cap`, is never overwritten. Fission yeast at L = 30 gets 859 nodes, `nmax` and
`mesh_reject` 42,858 and a cold cap of 21,429; Guzzo_noXA at L = 12 gets 481 / 963 / 801 nodes, 40,000 and 20,000.
The rule keeps the guess density at the L = 14 value, **57 nodes per unit length**; L = 7 runs at 114. Not scaled:
`identity.n_comparison_nodes`, `validity.n_dense`, `identity.relsup_nodes`. A mode-n shape uses only 1/n of the
nodes per piece; when (`guess_nodes` - 1) / n is below 50 the run prints a per-piece warning
(`SA:SEED:W_PIECE_NODES`). A bigger mesh limit lets failing cold calls run longer, and an attempt file at `nmax`
42,858 is about 19 MB.

## 7. Stops, reuse and resume

Stop reasons: `completed`, `budget_calls`, `budget_wall`, `budget_ring`, `disk_floor` (free space below
`checkpoint.min_free_gib`, default 8 GiB). A killed process can also be resumed.

The wrapper decides from the run folder (`sa_resume('status', ...)`, which writes nothing), in this order. A
`legacy` or `refused` folder stops the wrapper before the library step, so nothing is built for it:

| status | when | what the wrapper does |
|---|---|---|
| `new` | the folder is missing or empty | runs |
| `legacy` | a folder written by software 0.5 or earlier (no `internal/`) | refuses; the folder still loads, plots and traces into another folder (`opts.out`) |
| `refused` | the parameters, a scientific setting, the software, MATLAB, model, library or code version differ, the run stopped on an integrity error, or the folder is not empty and holds no run | refuses with every differing item, zero calls, the folder untouched: choose a new `run_dir`. When the parameters differ, the library loader's message (`SAP:LIBRARY` for fission yeast, `SAG:LIBRARY` for Guzzo_noXA, naming the parameter that does not match the library) may appear instead |
| `completed` | the run finished (also after compaction) | no solve and no replay; plot and traces only |
| `budget_exhausted` | stopped at `budget.max_calls` and the budget is not raised | no solve: raise `budget.max_calls` to continue |
| `resumable` | any other stop | resumes: every stored call is replayed and checked without solving, then the run continues |

Only `budget.*`, `checkpoint.*`, `report.*` and `import.source_dir` (a moved source) may differ on resume. Everything else, including the mode,
`lattice.step` and `lattice.points_per_axis`, is scientific. Edits to the plot, trace, wrapper, compaction, loader or
`post/` code never block a resume.

A finished run is checked like a resume. When package code that the run depends on (the `sa_*.m` files,
`seeded_atlas_run.m`, the model module, `solver/` and the private folders) has changed since the run started, even
a `completed` folder is `refused` (`E_CODE`), and the wrapper stops before the plot and the traces. To plot and
trace such a folder, call the tools directly; they make no solve, no replay and no status check, and write only
`plots/` and `traces/<name>/`:

```matlab
setup_atlas;                                   % from the package folder
seeded_atlas_plot(run_dir);                    % plots/roots_3d.fig, .csv, _key.mat and plots/sheets_2d/
seeded_atlas_trace(run_dir, struct('name', 'C1', 'root', 'k:1,4/R000001'));   % traces/C1/
```

Through the wrapper instead, add `settings_overrides.resume = true` and
`settings_overrides.resume_code_override = '<why the code changed>'` with the run's own settings and the same
MATLAB thread mode: the status is then `completed`, and still nothing is solved or replayed. A call cut off by a kill leaves an `interrupted` record that counts toward
`budget.max_calls`; the resume runs the same unit again.

**An import run and its source (since 1.0.0-rc2).** The source's identity is a stamp (SHA-256 over its
`settings_effective.mat`, `attempts.mat`, `queue_final.mat`, `journal.txt`, `compaction.mat` and the attempt-file
count). A resume of the follow-up recomputes it and refuses a changed source (`E_IMPORT`, after `E_LIBRARY` and
before `E_CODE`); each import call also checks the SHA-256 of its source attempt file (an integrity stop). A moved
source resumes with `import.source_dir` set to its new place. Resuming or compacting the source after the follow-up
has started changes the stamp and so blocks the follow-up's resume, by design; a `completed` follow-up keeps the
status `completed` (with a note), so it can still be plotted and traced. Plotting or tracing the source does not
change the stamp.

**A focus or import run from a killed source** (`import.source_stopped`) is identified by the SHA-256 of its
planned attempt files instead of the stamp, so resuming the source does not block it.

**The per-root pass and `root_pass_cap`.** With `root_pass.max_rounds` 10000 a `max_rounds` end is practically
unreachable. `root_pass.max_calls_per_round` is a user cost cap: when reached, the pass ends for the current search
pass (`root_pass_cap`); a later closure pass (when rescue added a node) may pose the rest. Imports do not change this.

**`from_run` across 1.0.0-rc2.** `sa_settings('from_run', ...)` keeps only the leaves you set, so a production run
made at rc2 from a pilot run made at rc1 resolves the new defaults (10000) unless you set the old ones.

Directly: `s.resume = true; s.budget.max_calls = 1000; seeded_atlas_run(p, manifest, s)`. For production after a
pilot, `[s, dropped] = sa_settings('from_run', pilot_dir, 'production')` keeps your own pilot settings and drops the
mode's (listed in `dropped`); set `s.run_dir` and the budgets.

## 8. Comparing with a reference

```matlab
ref = seeded_atlas_compare('from_myxo_sheet', '<myxo sheet>.mat', 0.15, 'myxo');   % or 'from_stageb_registry'
C = seeded_atlas_compare(run_dir, ref);   % writes run_dir/compare/myxo.csv and myxo_map.txt
```

Each node is classed agree (≤ 1e-6), band, different (≥ 1e-2), ours-only or reference-only. A reference never seeds
or steers a run; the classes are a report, not a commissioned tolerance.

## 9. What a result does not mean

A root is a steady state that passed the checks at these tolerances. It is not shown to be stable, two neighbouring
roots are not shown to lie on one branch, and the similarity tolerance (1e-3) is provisional. A node without a root
is not proof that none exists, and a call cut by the cold mesh cap is not evidence either. A trace's end state is not
shown to be an attractor. The 3-D plot's third coordinate is a display coordinate; it is not shown to separate
families.

## 10. Commissioned scope (at R2025b Update 4, maca64)

The pipeline is commissioned only at the settings it was validated at, and only on MATLAB R2025b Update 4 on macOS
arm64 (maca64); anything else (another L, other parameters, a new lattice, another MATLAB release or machine) is a new
acceptance test. The current list is in [MODELS](MODELS.md); the evidence (milestone reports and reviews) is in the
development repository, not in this package.

**What 0.12 was and was not tested on (milestone 8 of the development repository).** The new shapes and the floor rule are bit-identical to 0.11 at
the default `seed_modes` (fission yeast L = 7 box, Guzzo_noXA K2' box, explicit MPSS targets). The mode-n libraries
were built at fission yeast L = 30 (res 15, the L = 7 densities, modes 3 to 11 found pieces, mode 12 was empty) and for
Guzzo_noXA at L = 6, f = 1 (mode 3: 4 of 625 cells). Operating notes from that L = 30 work:

- **One pool worker gives no parallelism but doubles the MATLAB hours** (the client process and the worker both
  run): the L = 30 base build took 3,600 s per process, about 2.0 MATLAB-h, and each MPSS rescue pass about 13 min.
  Budget hours as processes times wall time, not wall time alone.
- **HPSS covers little at L = 30 / res 15**: 16 of 225 cells (the base libraries HSS 225/225 and MPSS 70/225). Families
  paired with the half-peak are therefore thin, and "no root" for them says little.
- **The scaled cold cap cut most mode-3/4 calls**: in the L = 30 coupled test the scaled cap (21,429) cut 74 of 81
  calls; five were accepted, four roots, none of mode 3 or 4. Absence of a mode-3/4 root there is a search result under
  the cap. **`cold.mesh_cap` is the lever** if those roots are wanted: a value you set wins over the derived one
  (`Inf` = no cap; it must not exceed `solver.mesh_reject`).
- Any run at L = 30 remains uncommissioned: the milestone tested the machinery, not validated results.

**What 1.0.0-rc2 was and was not tested on (milestone 10 of the development repository).** With the four cap leaves pinned at 8 / 4 / 16 / 100 and no
import, rc2 is bit-identical to rc1 on the fission-yeast L = 7 box and the Guzzo_noXA K2' box. Import was tested on
the fission-yeast L = 7 box (same grid and a 2/3-step grid), with kill and resume, and on 15 nodes of the L = 30 run
(re-solve time per mesh size); the mirror symmetry was checked on stored roots of both models. The L = 30 follow-up
run itself was not run; its cost is an extrapolation (milestone-10 packet section 7, development repository).

## 11. Operating notes

- **Generated solver copy.** The pipeline calls a traced copy of your MATLAB's own `bvp4c` that records every
  proposed mesh (`solver/sa_bvp4c_traced.m` and `solver/private/`, with the generation record
  `solver/sa_bvp4c_generated.txt`). These files are not in git and no MathWorks file is shipped. Make them with
  `setup_atlas` (from the package folder); `run_seeded_atlas.m` calls `setup_atlas('quiet')` before every run and
  remakes the copy when it is missing or does not match your MATLAB. `setup_atlas('check')` verifies only and prints
  one line, e.g. `setup_atlas: OK  R2025b (25.2.0.3150157) maca64  bvp4c 597839cef613  toolbox=reference
  deval=replica  selfcheck=pass(7/7)`, or `setup_atlas: FAIL <check>: <reason>`.
- **Other MATLAB releases (release policy).** Any release is accepted whose `bvp4c.m` the generator can patch and
  whose patched copy passes the self-check: on small fixed problems the copy must give exactly the same results as
  your untouched `bvp4c` (with a same-session control), its mesh hooks must fire where the pipeline expects, its text
  must equal your `bvp4c.m` except the renamed function line and two inserted lines, and the warning and error
  identifiers the pipeline relies on must exist. A toolbox equal to R2025b Update 4 is labelled `reference`, any
  other `self-checked`. A `.p` or MEX `bvp4c`, or one whose anchor lines moved (a rewrite), is refused with a named
  check (`E_TOOLBOX`, `E_ANCHOR`, `E_CLOSURE`, `E_SELFCHECK`); there is no mode without the mesh trace. The
  package's fast `deval` replica is used only on the reference `deval.m` (`deval=replica`); otherwise your own `deval`
  is used (`deval=toolbox`). The self-check certifies the patch, not the numerics: another release may find other
  roots, call counts or meshes. Tested only on R2025b Update 4, macOS arm64; other releases are self-checked, not
  tested; the minimum release is unverified.
- **MATLAB updates and other machines.** Every run records the MATLAB version (with the Update), the architecture
  and the toolbox hashes. A run cannot be resumed after a MATLAB Update or on another architecture (`E_VERSION`); a
  run started on another `bvp4c` is refused, and the run status says "the MATLAB toolbox changed since this run
  started". Finish a run on the MATLAB that started it, or start a new `run_dir`. Finished folders still load, plot,
  trace and compare.
- **Do not run the shipped tests while a run is in flight in the same package folder.** The tests that make traced
  copies (T61p, T102) work only in scratch folders and never touch the installed copy; the note is a precaution for
  the tests' scratch under `data/experiments/`.
- The uncoupled script `run_atlas_surface_pipeline.m` always reruns its own trace stage (`force_trace = true`); the
  wrapper never runs that stage.
