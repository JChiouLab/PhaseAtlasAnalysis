# Models, parameters and tested settings

The two models of the paper, in the paper's notation, with the names the code uses. Both are mass-conserving
reaction-diffusion models on a one-dimensional cell, x ∈ [0, L], with no flux at the two poles. The equations below
are written from the code (`src/models/guzzo_noxa/gnx_reaction.m`,
`src/solvers/atlas_surface_pipeline_dependencies/ode_sys_v7.m`) in the paper's symbols.

## The M. xanthus model (`guzzo_noxa`)

The model of Guzzo et al. (2018) reduced to RomR and MglB, with MglA eliminated through the effective activity
a_eff ("noXA"). The reference script of the reduction is kept byte for byte as
`src/pipelines/seeded_atlas_standalone/post/guzzo_noxa/Guzzo_noXA.m`; it is never executed.

```text
∂t rm = Dm ∂xx rm + F_R          ∂t rc = Dc ∂xx rc − F_R
∂t bm = Dm ∂xx bm + F_B          ∂t bc = Dc ∂xx bc − F_B
a_eff = rm / (1 + (bm/Ka)^2)
F_R   = (kr + krb bm) rc − dr rm
F_B   = (kb + kbb bm) bc − db bm / (1 + bm/K) − f · dba · a_eff · bm^2
```

Conserved: MR = ∫(rm + rc) dx and MB = ∫(bm + bc) dx. Concentrations are in units of the mean total densities
(`rho_R`, `rho_B`, both 1 in the paper); x in µm, t in min.

| code | paper | code value | Table 1 | meaning |
|---|---|---|---|---|
| `Rm Rc Bm Bc` | rm rc bm bc | | | membrane / cytosolic RomR and MglB |
| `kR` × `r_speed` | kr | 0.028624 × 0.35 = 0.0100 | 0.01 | basal RomR attachment |
| `kRB` × `r_speed` | krb | 1.0423 × 0.35 = 0.365 | 0.36 | RomR attachment promoted by membrane MglB |
| `dR` × `r_speed` | dr | 0.094962 × 0.35 = 0.0332 | 0.033 | RomR detachment |
| `kB` | kb | 8.1255 | 8.1 | basal MglB attachment |
| `kBB` | kbb | 3.5216 | 3.5 | cooperative MglB self-recruitment |
| `dB` | db | 5.4381 | 5.4 | MglB detachment |
| `K` | K | 0.11440 | 0.11 | saturation of the inhibition of MglB detachment |
| `dBA` × `dba_factor` | dba × f | 210.01 × f | 210 (f = 1) | MglB removal by a_eff |
| `KA` | Ka | 0.14714 | 0.15 | inhibition constant of a_eff by MglB |
| `DM` | Dm | 0.2 | 0.2 | membrane diffusion (µm² min⁻¹) |
| `DC` | Dc | 120 | 120 | cytosolic diffusion (µm² min⁻¹) |
| `L` | L | 6 | 6 | cell length (µm) |

Pipeline parameters (`par`): `L`, `r_speed` (multiplies `kR`, `kRB`, `dR`; 0.35 in the paper),
`dba_factor` (the signal factor f), `rho_R`, `rho_B`, `dn_half_range` (the half-cell library grid spans
±`dn_half_range` in `dn`), plus the library grid `res`. The paper's setting, with f = 0.02, 0.05, 0.08, 0.1 or 0.35
(Fig. 2E):

```matlab
par = struct('L', 6, 'r_speed', 0.35, 'dba_factor', 0.35, 'rho_R', 1, 'rho_B', 1, 'dn_half_range', 1.8); res = 25;
```

Plane coordinate: `dn = (n_left − n_right) / (L/2)`, the difference of the halves' mean densities; the paper's
fraction is Δn/M = `dn` / (2ρ). Both halves stay positive when |`dn`| < 2ρ. Uniform states do not depend on
`r_speed`; profiles with diffusion do.

## The S. pombe model (`original_polarity`)

```text
∂t u = Dm ∂xx u + (kgg u^2 + kg) v − dg u − kpg u p
∂t q = Dm ∂xx q + kpg u p − dgp q
∂t v = Dc ∂xx v − (kgg u^2 + kg) v + dg u + k* w
∂t w = Dc ∂xx w + dgp q − k* w
∂t p = Dc ∂xx p − kpg u p + dgp q
```

Conserved: MG = ∫(u + q + v + w) dx (GEF) and MP = ∫(p + q) dx (PAK); x in µm, t in s.

| code | paper | code value | Table 2 | meaning |
|---|---|---|---|---|
| `u q v w p` | u q v w p | | | membrane complexes u, q; cytosolic GEF v, phosphorylated GEF w, cytosolic PAK p |
| `a` | kgg | 2.38/4 = 0.595 | 0.59 | cooperative recruitment of GEF by the membrane complex |
| `b` | dg | 0.185 | 0.19 | detachment of the membrane complex |
| `c` | kpg | 151.6/2 = 75.8 | 76 | binding of Pak1 to the membrane complex |
| `d` | dgp | 0.552 | 0.55 | disassembly of the Pak1-bound complex |
| `e` | k* | 0.828 | 0.83 | dephosphorylation of cytosolic GEF |
| `f` | kg | 1 | 1 | basal attachment of cytosolic GEF |
| `Dm` | Dm | 0.05 | 0.05 | diffusion of u, q (µm² s⁻¹) |
| `Dc` | Dc | 5 | 5 | diffusion of v, w, p (µm² s⁻¹) |
| `M_BG` | MG | `2*L` | 2L | total GEF |
| `M_P` | MP | `2*0.246*L` | 2ρL | total PAK; the code uses ρ = 0.246 (Table 2: 0.25), and ρ = 0.246/3 for the Ashbya set |
| `L` | L | 7 (Ashbya 30) | 7 (Ashbya 30) | cell length (µm) |

The paper's settings, as in `run_seeded_atlas.m`:

```matlab
par = struct('a', 2.38/4, 'b', 0.185, 'c', 151.6/2, 'd', 0.552, 'e', 0.828, 'f', 1, ...
             'Dm', 0.05, 'Dc', 5, 'L', 7, 'M_BG', 2*7, 'M_P', 2*0.246*7);              % Fig. 3, S3, S4
par = struct('a', 2.38/4, 'b', 0.185, 'c', 151.6/2, 'd', 0.552, 'e', 0.828, 'f', 1, ...
             'Dm', 0.05, 'Dc', 5, 'L', 30, 'M_BG', 2*30, 'M_P', 2*(0.246/3)*30);       % Fig. 5, S5 (M_P = 4.92)
```

Write `par` with these expressions: a library is reused only when its stored parameters equal `par` bit for bit.
Plane coordinates: `d_G` = ΔGEF/MG and `d_P` = ΔPAK/MG, so the paper's ΔPAK/MP = `d_P`/ρ.

The half-cell library holds the three shapes: no-peak (`HSS`), half-peak at the pole (`HPSS`) and full-peak
(`MPSS`). The full-peak half is constructed from a half-peak solution on a quarter cell, reflected about its peak; it
is not solved directly as an interior-peak problem.

## Tested settings

All checks below were run on MATLAB R2025b Update 4, macOS arm64. Another release is accepted when the traced
`bvp4c` passes its self-check ([user guide](USER_GUIDE.md), section 11), but no result on it was tested.

- **Paper settings.** M. xanthus: L = 6, `r_speed` 0.35, library grid 25 × 25 (`dn_half_range` 1.8), grid step 0.15,
  `dba_factor` 0.02, 0.05, 0.08, 0.1 and 0.35 (production runs). S. pombe: L = 7 and L = 30 (Ashbya set) with the
  parameters above.
- **Other tested settings.** M. xanthus at the same L, `r_speed` and grid with `dba_factor` 0, 0.07 and 1: the
  computed states reproduce an independent reference computation where it has usable points, and runs restart
  correctly after the process is killed (`kill -9`) and resumed. S. pombe at L = 14: continuation with rescue and
  restart after a kill.
- **Full-model simulations** reproduce reference integrations: S. pombe bit for bit at L = 7.9 (ρ = 0.24621), and M. xanthus within
  1e-13 of the reference script `Guzzo_noXA.m` at 129 mesh points.
- **Not shown by any check:** that every quasi-steady state at a point was found, that neighbouring states lie on one
  branch, or their stability. A point without a solution is not proof that none exists.
- **Other L or parameters** run with the same code, but the settings (tolerances, mesh sizes, shape thresholds) were
  chosen at the settings above; start with a pilot run. Above L = 14 (S. pombe) or L = 6 (M. xanthus) the mesh and
  node counts are scaled with L (USER_GUIDE section 6, "The floor rule").
