# Seeded coupled atlas: operating contract (specification)

The product contract of the seeded coupled atlas, as adopted in the development repository (operating contract
version 1.1, 2026-10-12, with later milestone addenda), for software 1.0.0-rc2. The development history, decisions,
milestone packets, reports and reviews are not part of this package; where this page names them, they are in the
development repository. Usage is in the [user guide](USER_GUIDE.md); models and the commissioned scope are in
[MODELS](MODELS.md).

## Deliverable and claim ceiling

Deliver a MATLAB pipeline for bounded discovery of valid coupled/forced roots,
with smoke/pilot/production (and focus) modes, inspectable profiles and diagnostics, shared
root identity, saved effective settings and restart. The user can inspect a
completed pilot and tune a new invocation; no AI or manual intervention is
required between stages within an invocation.

The coupled pipeline reads seed libraries written by the separate uncoupled
atlas pipeline, which is expanded to run Guzzo_noXA as well as the original
model. Each model sits behind a model module (equations, fields, pools,
normalization, library reading, model-specific validity details). Parameters are
user-chosen; each model–pipeline pairing stays `experimental` until its own
commissioning checks pass, and no tolerance or taxonomy is inherited across
models without a check.

Formal closure, completeness, independently confirmed multiplicity, stability,
branch connectivity and PDE correspondence are separate claims, not first-product
acceptance requirements. Rendering and PDE integration remain later work.
Unattended acceptance requires a fresh invocation; historical replay alone is
insufficient. Stage-B closure is not a prerequisite for this product.

## Numerical and root acceptance

The original-model adapter must apply one acceptance boundary before roots enter
the accepted store or become eligible parents: finite values, defining-equation
residuals, boundary/seam conditions, mass constraints and positivity within
prospectively declared tolerances. Preserve source/flux consistency diagnostics
and specify their applicable acceptance tests and near-zero scaling. A low solver
residual alone is insufficient; clamping a reaction evaluation cannot establish
admissibility of the original profile. Invalid returns remain diagnostic records,
not accepted parents. Recovery and exports use the same saved validity semantics;
missing validity evidence fails closed rather than relying on plotting filters.

Set tolerances and representative refinement checks for the actual model,
normalization and formulation in the implementation/run packet. Do not import
myxo thresholds or infer continuum accuracy from solver convergence.

All discovery routes share one root store. Compare all physical fields with a
declared scale-aware distance at the same target; retain distances, profiles and
near-threshold cases. Family names describe shape/provenance and cannot create
or remove a root. Preserve distinct valid returns and flag large profile changes
without claiming branch continuity. Reuse the audited M1 kernel only within its
validated software contract; its scientific similarity tolerance is not yet
commissioned.

## Search and performance

Retain library-seeded cold search, continuation, rescue and exhaustive cold search
as a reference option. Since 1.0.0-rc2 (milestone 10 of the development repository) a declared `import` route may
seed a run from the roots of a finished run folder of the same model and
parameters: each root (and, by default, its mirror image) is re-solved once at the
nearest node, under the run's own validity boundary, with the source identified by
a stamp; an import is evidence of its own call only, and the source run is never a
hidden seed (the route, stamp and source path are recorded and reported). Since
milestone 11 (also 1.0.0-rc2) a fourth mode, `focus`, imports a declared selection of those roots
(by id or observed family label) as the only seeds and fences the search to the
selected families: roots outside the fence are kept and reported but never
stepped from; a killed source is accepted only on an explicit setting, with
its calls lacking an end line ignored. Since milestone 12 of the development repository (also 1.0.0-rc2; owner
ruling 2026-11-10) a `focus` run uses no uncoupled library: the mode itself is the
explicit configuration; its identity field scale and family lookup are taken
from the source run's stored effective settings (recorded with the source
identity, checked on resume), and no census, cold targets or ring are
computed. Add an optional anchor-first strategy using a small,
declared set of direct anchors or library seeds. Direct anchors supplement the
library route; library-free use must be explicitly configured, never silently
triggered by a missing or malformed library. Report library seedability separately
from anchor/search reachability. Large parent-free ensembles remain deferred.
The plan's central-nine option remains a specific library-seed selection policy.

Integrate bounded neighbor rescue, including appropriate single-parent retry and
pair-based routes, with the shared attempt history. Revisit newly eligible work
within the configured budget; a failed march must not silently label the rest of
a row refused or nonexistent. Distinguish invalid, refused, unattempted and
policy-skipped targets. Queue exhaustion certifies only the configured search.

Parallel execution is optional after the serial path works. Use independent
work units where dependencies permit, one durable writer and deterministic root
registration. Verify serial/parallel agreement on representative cases before
relying on the optimization; document any permitted numerical differences.

Compare strategies using valid distinct roots, valid target coverage, missed
reference roots, total time/calls and refusal reasons. Use a declared target set
and equivalent resource accounting, including failed solves and rescue. Historical
coupled roots may be comparison references after fresh discovery, not hidden
seeds. Extra coverage is not an all-roots claim.

## Seed shapes, mode libraries and point-count scaling (since software 0.12)

Contract only; usage is in the [user guide](USER_GUIDE.md); the design packet and the review evidence (R17) are in
the development repository.

- **`seed_modes`.** A shared scientific setting: a nonempty cell row of left-drawn half shapes (`0`, or `n+` / `n-`,
  n >= 1) or `'all'`. Each listed shape seeds both halves; the right half's drawn name is the mirror of the shape (sign
  flips for odd n). Families are the ordered pairs, filtered by the census and `families`; there is no pairing
  restriction. Default `{'0','1+','2-'}`: at the default the 0.11 code path runs and results are bit-identical
  (families, census, targets, order, rows, stores). `'all'` resolves from the library folder to the shapes it offers and
  stores the resolved list; each skipped `'empty'` mode file is named in a printed line. `1-` (HPSS guess reversed) and
  `2+` (MPSS stitch rotated by a quarter) are made on load from existing files; `n+` / `n-` (n >= 3) need
  `branch_MODE<n>_latest.mat`. Shared code names no half type: the shape-to-type map, type order and mirror rule are
  model-module members (`seed_shape`).
- **Refusals, at zero calls.** Unknown or repeated name, `'all'` mixed with names (`SA:SEED:E_SHAPE`); a needed mode file
  missing (`SA:SEED:E_FILE`, naming the file and the `max_mode` that builds it); a named mode file with `run_status`
  `'empty'` (`SA:SEED:E_EMPTY`); a library file at other parameters (`SAP:LIBRARY`) or on another grid.
- **Pre-run table.** Printed before the first call (zero solves): per seedable family the expected drawn label, seeds,
  declared cold calls and rough hours at the model's measured cold rate (fission yeast 3.7 s, Guzzo_noXA 3.0 s); pairs with
  no seed listed as `0 (no cell with both halves)`; families with a new shape tagged as having no measured rate. The hours
  count cold calls only and are an order of magnitude, not a budget.
- **Mode libraries.** Built by `atlas_build_libraries(model, par, res, out_dir, struct('max_mode', n, 'pool_workers', k))`
  after the three base libraries, with an explicit nonempty `out_dir` (required when `max_mode` > 2): a half-peak piece on
  L_half / n at 1/n of the half masses, searched by the existing sweeper unchanged, stitched n times, accepted by the
  legacy gates plus a trough-end slope gate; the grid is copied from the sibling files and `par` is the original. The
  n = 2 branch reproduces the stored MPSS libraries bit for bit. The build stops after the first empty mode (file still
  written, `run_status` `'empty'`): a cost rule; empty is a search result, never non-existence. A mode file no active
  shape needs is not loaded, does not change the library fingerprint and does not block a resume.
- **Floor rule.** `scaling.reference_length` is a model-layer constant (fission yeast 14, Guzzo_noXA 6; in the model
  digest; package default Inf = off). At L <= reference length nothing changes. Above it, for each of
  `model_options.guess_nodes`, `continuation.guess_nodes`, `model_options.shape_nodes` (1 + ceil((N - 1) L / L_ref)) and
  `solver.nmax`, `solver.mesh_reject`, `cold.mesh_cap` (ceil(N L / L_ref)), whose source is `package` or `model`, the value
  is scaled and its source becomes `derived`; user values win; scaling that leaves mesh_reject > nmax or a cold cap above
  mesh_reject is refused. Resulting guess density 114 nodes per unit length at L = 7 and 57 at L >= 14. Per-piece
  warning `SA:SEED:W_PIECE_NODES` when (`guess_nodes` - 1) / n < 50.
- **`targets.policy = 'per_seed_centre'`.** per_seed's targets with the cold calls ordered nearest d = 0 first; a run
  option, never a default.
- **Claim ceiling.** Bit identity at defaults and the L = 30 library and coupled sizing are tested
  (R17: PASS WITH REPAIRS). A mode-n library that is empty, or a mode-3/4 coupled root not found under the scaled cold
  cap (L = 30: 74 of 81 calls cut at 21,429; `cold.mesh_cap` is the lever), is a search result. HPSS at L = 30 / res 15
  covers 16 of 225 cells. One pool worker doubles MATLAB-h without parallelism. L = 30 stays uncommissioned.

## MATLAB releases (since 1.0.0-rc1)

- The traced `bvp4c` copy is generated from the user's own MATLAB and accepted only when it passes the self-check
  (bit-identical results to the untouched `bvp4c` on fixed problems, mesh hooks where the per-call consistency rule
  expects them, a text equal to the toolbox file except three lines, the identifiers the pipeline relies on). The
  R2025b Update 4 hashes are a reference record (`toolbox=reference`); any other toolbox is `self-checked`. An
  unpatchable `bvp4c` is refused with a named check; there is no trace-less mode.
- The fast `deval` replica runs only on the reference `deval.m` and `private/ntrp3h.m`; otherwise the toolbox `deval`
  is used.
- Every run records the software version, the MATLAB version (with the Update), the architecture and the thread mode
  (never waived on resume) and the toolbox hashes in the model version.
- Claim ceiling: the self-check certifies the patch, not a release's numerics. The commissioned scope holds at R2025b
  Update 4, maca64 only; other releases are self-checked, not tested; the minimum release is unverified.

## Sources

The adopted direction and every milestone's evidence (packets, reports, independent reviews) are recorded in the
development repository (private). This page is maintained there and copied into each export.
