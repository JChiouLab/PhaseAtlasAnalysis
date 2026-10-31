# Changelog

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
