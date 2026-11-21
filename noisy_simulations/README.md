# Noisy simulations (manuscript Fig. 4)

Stochastic simulations of the fission-yeast polarity model (five species u, q, v, w, p; u and q on the membrane),
with Ornstein-Uhlenbeck noise on the membrane species. This folder contains the code that produced the 50 noisy
runs, the noise-free control and the long run of Fig. 4. It is Python; the rest of this repository is MATLAB and is not needed here.

## Credit

The solver and noise model, `gef_pak_noise/`, are the package of C. Hinds, used here with minor edits to comments
and docstrings only (the code is unchanged). Its equation references are to Dao et al. (2026). `run_comparison.py`
is an example script of that package with placeholder parameters; it needs matplotlib and is not used for Fig. 4.

## Files

| file | role |
|---|---|
| `gef_pak_noise/` | the noise package (C. Hinds): model (`model.py`), grid and Laplacian (`grid.py`), OU noise (`noise.py`), Crank-Nicolson solver with fixed-point iteration (`solver.py`) |
| `noise_driver.py` | one run per call: the manuscript parameters, initial condition, seeds, the mass-conserving noise (`MassFixSolver`) and the observables |
| `run_fig4_batch.sh` | the batch of 50 replicates at A = 0.001 and the A = 0 control, then the classification |
| `analyze_runs.py` | labels each run's time course as oscillation or monopolar epochs and classifies the runs (criterion in its docstring) |
| `tests/test_reproduce.py` | reruns the first 40 s of two Fig. 4 runs and compares them with the stored production output (`tests/reference_first40s.npz`) |

## How the Fig. 4 panels were produced

- **Fig. 4A, B:** the 50 replicates at A = 0.001 (rep 0-49) with the defaults (T = 4000 s), from
  `run_fig4_batch.sh`. The plotting code is not included.
- **Fig. 4C:** one long run, rep 100, `python3 noise_driver.py --out runs --A 0.001 --rep 100 --T 14400`
  (240 min). It is classified with a 5-min monopolar threshold:
  `FIG4_T_LONG=300 python3 analyze_runs.py <runs folder> <output prefix>`. Without `FIG4_T_LONG` the threshold is
  3 times the noise-free half-period (141 s), i.e. 423 s.
- The state is rho = M_P/M_BG = 3.444/14 = 0.246 (Table 2 rounds this to 0.25).

## Requirements

Python 3 with numpy and scipy. Fig. 4 was produced with Python 3.11.8, numpy 1.26.4 and scipy 1.12.0
(`requirements.txt`); with these versions the test reproduces the stored output exactly.

## Quick start

```sh
python3 tests/test_reproduce.py                      # about 10 s; prints PASS
python3 noise_driver.py --out runs --A 0.001 --rep 0 # one Fig. 4 replicate, about 20-30 min
sh run_fig4_batch.sh runs 4                          # all 51 runs, 4 at a time, then the classification
```

Each run writes `A<A>_rep<rep>.npz` with the observables every 0.5 s (`t`, `dG`, `dP`, `dn2_u`, `uL_max`,
`uR_max`, total masses `M_GEF`, `M_PAK`, `min_conc`), u profiles every 2 s (`snap_t`, `snap_u`), the final full
state and a `meta` record (parameters, seeds, wall time, status).

## Model and settings

- Kinetics a = 0.595, b = 0.185, c = 75.8, d = 0.552, e = 0.828, f = 1; diffusion D_m = 0.05 (u, q),
  D_c = 5 (v, w, p); cell length L = 7 with no-flux ends; total GEF 14 and total PAK 3.444 (rho = 0.246).
- Grid 200 nodes (dx = 7/199), Crank-Nicolson with fixed-point iteration, dt = 0.0025 s, T = 4000 s.
  At dt = 0.01 and 0.005 the fixed-point iteration fails when the first peak forms.
- Noise: independent OU fields on u and q at every node, correlation time τ_c = 1, strength A = σ²τ_c = 0.001,
  advanced by Euler-Maruyama and averaged over each time step.
- **Mass conservation (`--mode massfix`, the default and the Fig. 4 setting).** The spatial mean of each noise
  increment is removed before it is added, so total GEF and PAK are conserved exactly in every realization
  (to about 1e-10). In the package as written (the uncorrected mode, `--mode uncorrected`) the noise has zero mean
  only in expectation, and the totals drift as a random walk.
- Initial condition: the homogeneous steady state plus a 1 % random perturbation (seed = rep).
  Noise seed = 100000 + 100·round(10⁴A) + rep.

## Observables

- dG, dP: (left − right) half-cell amounts of GEF = u + q + v + w and PAK = p + q, divided by the total GEF (14).
  Normalised to each protein's own total instead: ΔGEF = dG and ΔPAK = dP/ρ with ρ = 0.246
  (normalised to half of each total: twice these).
- dn2_u = ∫₀ᴸ cos(2πx/L) u dx, the shape axis of Figs. 3 and 4.
- uL_max, uR_max: the maximum of u in each half (the manuscript plots u divided by the mean GEF density 2).
