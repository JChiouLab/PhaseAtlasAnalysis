# Changelog

## Unreleased (2026-11-18)

- **Added `noisy_simulations/`** (Python): the code of the manuscript's Fig. 4 noisy simulations. It contains Clara
  Hinds' solver and noise package `gef_pak_noise/` (unchanged), the run driver `noise_driver.py` with the
  mass-conserving noise used for Fig. 4, the batch script, the run classifier and a reproduction test that matches the
  stored production output bit for bit. Added directly to this repository, not exported from the development
  repository; its files are listed in `EXPORT_MANIFEST.txt` with origin `direct-add`.

## 1.0.0-rc2 (2026-11-10)

New features (development milestones 10 and 11; details in `docs/USER_GUIDE.md` section 3):

- **Roots of an earlier run as seeds (import).** Set `import_from` in `run_seeded_atlas.m` to a finished run folder of
  the same model and parameters. Each of its roots, and by default its left/right mirror image, is solved once at the
  nearest node of the new grid (the grid may be the same, finer, or not nested) and then serves as a starting root
  like a cold root. Optional: skip steps the earlier run already tried and refused (`import.source_edges`
  `'same_or_stronger'` or `'all'`; same grid, solver settings and model version only; default `'none'`). The earlier
  run is identified by a stamp, and a resume refuses a changed source (`E_IMPORT`).
- **Focus mode** (`run_mode = 'focus'`). Pick roots of an earlier run, by id (`import.root_ids`, a list or a CSV with a
  `root_id` column) or by family label (`import.families`); only those are seeded (no cold calls), and the search is
  fenced to their families and the mirror families: roots of other families are kept and reported but never stepped
  from. A budget and `import_from` are required.
- **The family fence** (`continuation.families`) can also be set in the other modes (off by default there).
- **Killed sources.** A run that was killed can serve as a source with `import.source_stopped = true`; its unfinished
  calls are ignored, and it is identified by its attempt files, so resuming it later does not block the new run.
- **Source-edge reuse requires the same model version** (the corrector and solver code), since a refusal depends on it.
- **`max_mode` in the edit block of `run_seeded_atlas.m`** builds mode-n seed libraries beside the base files; with
  `seed_modes = 'all'` and no mode files present, the run warns (`SA:SEED:W_NO_MODE_FILES`).
- **Compaction of a cleanly stopped run** (stopped on disk, budget or wall time, not killed):
  `seeded_atlas_compact(run_dir, struct('stopped_ok', true))`. The run can no longer be resumed; it still plots,
  traces and serves as an import source.

Default changes:

- `root_pass.max_rounds` 8 -> 10000, `multiseed.max_roots_per_node` 4 -> 10000, `continuation.max_lineages`
  16 -> 10000, `continuation.max_calls_per_lineage` 100 -> 10000: the per-root pass runs until a round adds nothing,
  and lineages are no longer cut by these caps. A pilot's projection prints the cap column as "unbounded".

Compatibility:

- With those four settings pinned at 8 / 4 / 16 / 100 and no import, rc2 gives bit for bit the same results as rc1
  (checked on the fission-yeast L = 7 box and the Guzzo_noXA K2' box). At the new defaults, runs differ wherever an
  old cap bound.
- `model.version` is unchanged. Run folders record `seeded_atlas_standalone/1.0.0-rc2`; folders written by rc1 still
  load, plot, trace and serve as import sources, but cannot be resumed by rc2.

### Update 2026-11-11: focus runs without an uncoupled library (development milestone 12)

The version number stays 1.0.0-rc2. Details in `docs/USER_GUIDE.md` section 3 ("Focus mode") and section 7.

- **Focus needs no library.** A `run_mode = 'focus'` run neither builds nor reads an uncoupled library: the selected
  roots of the earlier run are its only seeds. `library_dir`, `res`, `max_mode` and `seed_modes` are ignored (the
  wrapper prints `library none (mode focus)`). The two things a run used to take from the library come from the source
  run instead: the identity field scale and the family lookup (mode name to library family; the source's
  `seed_modes` replaces the run's own). No seedability census and no cold targets are computed; `routes.cold = true` in
  focus is refused (`SA:IMPORT:E_FOCUS_COLD`).
- **No ring in focus.** Without a library there is no library range, so every lattice node counts as interior:
  `search.ring_last`, `rescue.ring` and `budget.ring_calls` have no effect in focus, and the outer band is searched
  and rescued like the rest.
- **Field-scale settings must match the source.** The settings the field-scale rule reads (`lattice.anchor`,
  `identity.field_scale_floor`, `model_options.hss_branch`; declared by a model module in the new optional member
  `field_scale_settings`) must equal the source's; a difference, or a source that lacks the stored items (written
  before software 0.11), is refused before any folder is made (`SA:IMPORT:E_FOCUS_SOURCE`, naming the setting or the
  missing items).
- **Recorded and checked on resume.** The run stores where its field scale and family lookup came from (source path,
  SHA-256 of the source's `settings_effective.mat`, the import identity); a resume derives them again and refuses on
  any difference or an unreadable source (`SA:RESUME:E_FOCUS`). A focus folder started with a library by the first
  1.0.0-rc2 build cannot be resumed by this one.
- Other modes (smoke, pilot, production, import without focus) are unchanged.


## 1.0.0-rc1 (2026-10-31)

The first export of the two pipelines as one standalone package, from the development repository's software 0.12
(seeded coupled atlas) and its uncoupled atlas pipeline, with:

- **Release self-check of the traced `bvp4c`.** The copy is made from the user's own MATLAB (`setup_atlas`) and
  accepted when it passes a self-check: bit-identical results to the untouched `bvp4c` on small fixed problems, mesh
  hooks where the pipeline expects them, a three-line text diff, and the warning and error identifiers the pipeline
  relies on. The former lock to one MATLAB release (R2025b Update 4) is gone; those hashes are now a reference record
  (`toolbox=reference`, otherwise `self-checked`). A generation record (`solver/sa_bvp4c_generated.txt`, ignored by
  version control) names the release, the toolbox hashes, the anchors and the self-check results.
- **`deval` replica only on the reference `deval.m`**; on any other release the user's own `deval` is used.
- **`provenance.arch`**: a run records the machine architecture; no resume across architectures (or MATLAB Updates).
- **Toolbox-change message**: a run started on another `bvp4c` is refused on resume, and the run status says "the
  MATLAB toolbox changed since this run started (bvp4c <old> -> <new>)".
- **`setup_atlas`** (path, generation, self-check; `setup_atlas('check')` prints a one-line verdict), and
  `run_seeded_atlas.m` remaking a missing or stale copy before every run.
- Run folders record `seeded_atlas_standalone/1.0.0-rc1`; folders written by 0.12 or earlier still load, plot and
  trace but cannot be resumed.

Earlier history: the development repository (private).
