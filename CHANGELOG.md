# Changelog

## 1.0.0-rc3

- **2-D sheets as PDFs in the style of Fig. 5.** `seeded_atlas_plot` (and so every run) now draws
  `plots/sheets_2d/<sheet>.pdf` and `all_families.pdf` instead of `.fig` files. Each sheet shows the grid points with
  a quasi-steady state, unit flux arrows from each state's sources, the nullclines and an example profile as an
  inset. The drawing is done by Python (`post/sa_plot_sheets.py`; numpy, matplotlib, contourpy) from
  `sheets_2d_roots.csv` and `sheets_2d.json`, which are always written. Without Python the command that draws the
  PDFs is printed; choose the Python with `opts.python` or `SA_PYTHON`. The 3-D `roots_3d.fig` is unchanged. New
  test `sa_t87_sheet_pdfs`.
- **`noisy_simulations/`** (Python): the code of the Fig. 4 noisy simulations: the noise package `gef_pak_noise/`,
  the run driver `noise_driver.py` with the mass-conserving noise, the batch script, the run classifier and a test
  that reproduces the first 40 s of two stored runs exactly.
- Documentation rewritten for readers of the paper: worked examples of the paper's three atlases, a figure map and
  a "Paper ↔ code" appendix in `README.md`; the models in the paper's notation in `docs/MODELS.md`.
- Licence: CC BY-NC 4.0 (`LICENSE`). Credits in `README.md` confirmed by all authors.
- Run folders written by 1.0.0-rc2 or earlier cannot be resumed with this version (the software version and code
  hashes differ); they can still be plotted, traced and used as `import_from` sources.

## 1.0.0-rc2

- **Solutions of an earlier run as starting points** (`import_from` in `run_seeded_atlas.m`): each quasi-steady state
  of a finished run of the same model and parameters, and by default its mirror image, is solved once at the nearest
  point of the new grid. Optionally, steps the earlier run already tried and refused are skipped
  (`import.source_edges`).
- **Focus mode** (`run_mode = 'focus'`): uses only chosen states of an earlier run (by id or by sheet label), needs no
  half-cell library, and restricts the search to their sheets (`continuation.families`).
- **Killed runs as sources** (`import.source_stopped = true`).
- **Higher-mode libraries** from the edit block (`max_mode`).
- **Compaction of a cleanly stopped run**: `seeded_atlas_compact(run_dir, struct('stopped_ok', true))`.
- **Default change:** `root_pass.max_rounds`, `multiseed.max_roots_per_node`, `continuation.max_lineages` and
  `continuation.max_calls_per_lineage` now default to 10000 (effectively unbounded), so continuation is no longer cut
  by these caps. Set the former values (8, 4, 16, 100) to reproduce an rc1 run.
- Run folders written by rc1 still load, plot, trace and serve as import sources, but cannot be resumed.

## 1.0.0-rc1

- First release of the two pipelines as one package.
- `setup_atlas` makes the traced `bvp4c` copy from your own MATLAB and accepts it after a self-check, so releases
  other than R2025b can be used; `setup_atlas('check')` prints a one-line verdict.
- A run records the machine architecture and the MATLAB release and refuses to resume on another or after a change
  of `bvp4c`.
