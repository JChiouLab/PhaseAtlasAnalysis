# Models and commissioned scope

This page lists the two models the package runs, their parameters, and where each model–pipeline pairing is
commissioned. It is condensed from the development repository's portability page and model pages; the evidence
(milestone reports, independent reviews) stays there.

## How a pairing is classified

A model–pipeline pairing is:

- `commissioned`: its declared checks passed at the named settings;
- `experimental`: under investigation, or used outside its commissioned settings;
- `unsupported`: a needed capability is absent.

Absence of evidence is not proof of impossibility, and nothing is inherited across models, parameter sets or
lengths without a check: a new L, new parameters or a new lattice is a new acceptance test. Every claim below holds
at **MATLAB R2025b Update 4 on macOS arm64 (maca64)** only. Another release is accepted by the package when the
traced `bvp4c` passes its self-check ([user guide](USER_GUIDE.md), section 11), but no result on it is commissioned.

## Commissioned scope (R2025b Update 4, maca64)

| model × pipeline | status and settings | not claimed |
|---|---|---|
| Guzzo_noXA × uncoupled atlas (library stages) | `commissioned` at the owner's reference setting only: L = 6 µm, `r_speed` 0.35, the 25 × 25 grid (`dn_half_range` 1.8) | other L, speeds or grids |
| Guzzo_noXA × coupled pipeline, cold path | `commissioned` at L = 6, s = 0.35, the 25 × 25 grid, `dba_factor` ∈ {0, 0.02, 0.05, 0.07, 0.35, 1}: accepted states are valid and reproduce the myxo atlas where it has usable cells | coverage, uniqueness, stability |
| Guzzo_noXA × coupled pipeline, continuation and restart | `commissioned` at the same setting, step 0.15, N = 241 (all six sheets) and N = 482 (`dba_factor` 0.35 production; `dba_factor` 1 `all_seedable` production, all reference roots re-found); restart shown by real `kill -9` runs | completeness, uniqueness, stability, connectivity; rescue reach (fails closed); the outer ring (\|Δn\| = 1.95) |
| Fission yeast × coupled pipeline, continuation with rescue | `commissioned` at L = 14 µm, the reference parameters, N = 401, box [3 3]; restart shown with real kills at L = 14 | other L (L = 7 is in use but not commissioned) |
| Both × coupled pipeline, modes and pilot reporting | `commissioned` at the tested settings: Guzzo_noXA L = 6, s = 0.35, step 0.15 (smoke, pilots, production); fission yeast L = 14 (smoke, pilot) | the production projection is an estimate |
| Both × coupled pipeline, route defaults (per-root pass, secant rescue, cold mesh cap 10,000, `per_seed` cold targets) | `commissioned` at fission yeast L = 14 (smoke, pilot, box production) and Guzzo_noXA (dBA × 0.35 production, N = 482; dBA × 1 production) | other L or parameters; under `per_seed` multiplicity is tested only by the later phases |
| Both × coupled pipeline, mesh trace and `sa_deval` replica | bit-identical to the earlier software at the tested settings (fission yeast L = 14 cold, refinement, continuation; Guzzo_noXA f = 1 production, K2′ box) | bit-identity across MATLAB releases or updates |
| Both × coupled pipeline, higher modes (`seed_modes`, mode-n libraries, floor rule; since 0.12) | at the default `seed_modes` and at or below the reference length (fission yeast 14 µm, Guzzo_noXA 6 µm) bit-identical to 0.11; the new shapes are exercised, not commissioned | L = 30 (uncommissioned) |
| Both × 1.0.0-rc1 (this export) | the release self-check accepts R2025b U4 (`reference`, `deval=replica`); the export is checked against 0.12 on two coupled smoke runs and the Guzzo_noXA library (export report) | any other release |

Plots and full-PDE traces are displays and integrations, not validated results; fission-yeast traces reproduced a
reference trace (Fig 3 trace C) bit for bit at the h04 setting, and Guzzo_noXA traces agree with Wen-Han's `Guzzo_noXA.m` to 1e-13 at
nx 129.

## Fission yeast (`original_polarity`)

Five fields (`u q v w p`) on each half of a cell of length L, two conserved pools (GEF and PAK). Parameters
(`par`): `a b c d e f Dm Dc L M_BG M_P` (rates, diffusivities, cell length in µm, total GEF and PAK masses). The
reference set in `run_seeded_atlas.m` and `run_atlas_surface_pipeline.m`:

```matlab
par = struct('a', 2.38/4, 'b', 0.185, 'c', 151.6/2, 'd', 0.552, 'e', 0.828, 'f', 1, ...
             'Dm', 0.05, 'Dc', 5, 'L', 7, 'M_BG', 2*7, 'M_P', 2*0.246*7);
```

Write `par` with these expressions: a library is reused only when its stored parameters equal `par` bit for bit.
The uncoupled seed library holds three half-cell families: HSS (flat), HPSS (half peak at the pole) and MPSS
(middle peak, constructed from a reflected quarter-domain half peak; never a direct interior-peak solve).

## Guzzo_noXA (`guzzo_noxa`)

A four-field myxobacterial polarity model (Wen-Han's reduction of Guzzo's model; a byte copy of his script ships as
`src/pipelines/seeded_atlas_standalone/post/guzzo_noxa/Guzzo_noXA.m`, used only as part of a trace key).

| index | field | meaning | diffusion |
|---|---|---|---|
| 1 | `Rm` | membrane RomR | `DM` |
| 2 | `Rc` | cytosolic RomR | `DC` |
| 3 | `Bm` | membrane MglB | `DM` |
| 4 | `Bc` | cytosolic MglB | `DC` |

```text
dRm/dt = DM Rm_xx + sR        dRc/dt = DC Rc_xx - sR
dBm/dt = DM Bm_xx + sB        dBc/dt = DC Bc_xx - sB
Aeff = Rm / (1 + (Bm/KA)^2)
sR   = (kR + kRB*Bm)*Rc - dR*Rm
sB   = (kB + kBB*Bm)*Bc - dB*Bm/(1 + Bm/K) - dBA*Aeff*Bm^2
```

Zero flux at both ends; x in µm, t in min. Conserved pools `N_R = ∫(Rm + Rc) dx`, `N_B = ∫(Bm + Bc) dx`, with mean
densities `rho_R`, `rho_B` (reference 1).

| name | value | meaning |
|---|---|---|
| `kR` | 0.028624 | basal RomR attachment |
| `kRB` | 1.0423 | RomR attachment promoted by Bm |
| `kB` | 8.1255 | basal MglB attachment |
| `kBB` | 3.5216 | cooperative MglB self-recruitment |
| `dR` | 0.094962 | RomR detachment |
| `dB` | 5.4381 | maximal saturating MglB detachment |
| `dBA` | 210.01 | MglB removal by Aeff Bm² |
| `K` | 0.11440 | saturation of MglB detachment |
| `KA` | 0.14714 | half-saturation of Aeff inhibition by Bm |
| `DM` | 0.2 | membrane diffusivity (µm² min⁻¹) |
| `DC` | 120.0 | cytosolic diffusivity (µm² min⁻¹) |

Pipeline parameters (`par`): `L` (µm), `r_speed` (multiplies `kR`, `kRB`, `dR`), `dba_factor` (multiplies `dBA`),
`rho_R`, `rho_B`, `dn_half_range` (the half-domain partition grid spans ±`dn_half_range`), plus the library grid
`res`. The reference setting is

```matlab
par = struct('L', 6, 'r_speed', 0.35, 'dba_factor', 1, 'rho_R', 1, 'rho_B', 1, 'dn_half_range', 1.8); res = 25;
```

Half-domain partition coordinate: `dn = (n_left − n_right) / Lh`, the difference of the halves' mean densities
(`Lh = L/2`); both halves stay positive when `|dn| < 2 rho`. Uniform states do not depend on `r_speed`; steady
problems with diffusion do.
