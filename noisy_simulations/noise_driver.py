"""Noisy simulations of the fission-yeast polarity model (manuscript Fig. 4).

One call runs one replicate (noise amplitude A, replicate index rep) of the five-species model
(u, q, v, w, p) with Ornstein-Uhlenbeck noise on the two membrane species u and q, and writes
<out>/A<A>_rep<rep>.npz. An existing output file is skipped unless --force is given.

The solver and the noise model are the package `gef_pak_noise` (in this folder; see README.md).
This driver adds:
  - the parameter set of the manuscript (L = 7; below), the grid, the initial condition and the seeds;
  - `MassFixSolver` (mode 'massfix', the default and the mode used for Fig. 4): the same noise fields,
    but the spatial mean of each noise increment is removed before it is added, so total GEF and total
    PAK are conserved exactly in every realization. The uncorrected mode (`--mode uncorrected`) runs the
    package solver exactly as written; there the noise is zero-mean only in expectation and the total
    masses drift;
  - the observables used in the manuscript (half-mass differences dG, dP, the shape coordinate dn2_u,
    the maximum of u in each half).

Defaults are the values used for Fig. 4: mode massfix, T = 4000 s, dt = 0.0025 s, A = 0.001
(the control run uses A = 0), output every 0.5 s, u snapshots every 2 s.

Grid. The package's Neumann grid puts N nodes at x_i = i * L_grid / N with mirror ghosts, so the no-flux domain
spans L_grid * (N - 1) / N. The driver passes L_grid = L * N / (N - 1), so the span is exactly L
(node-based, dx = L / (N - 1)).

Initial condition. The unique homogeneous steady state (HSS) plus a 1 % relative Gaussian perturbation
(seed = rep, so every A uses the same initial condition for a given rep), with the perturbation's
trapezoid-weighted mean removed per species so that the total masses are exact.

Noise seed = 100000 + 100 * round(1e4 * A) + rep.
"""
from __future__ import annotations

import argparse
import os
import sys
import time

import numpy as np

# make the package in this folder importable however the script is started (also under python -I)
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

from gef_pak_noise.grid import Grid                                                     # noqa: E402
from gef_pak_noise.model import ReactionParams, solve_homogeneous_equilibrium, reaction  # noqa: E402
from gef_pak_noise.solver import MCRDSolver, SimulationResult                          # noqa: E402
from gef_pak_noise.noise import OUNoiseField                                            # noqa: E402

# ---------------------------------------------------------------- parameters (manuscript, L = 7)
L = 7.0
PARAMS = ReactionParams(a=0.595, b=0.185, c=75.8, d=0.552, e=0.828, f=1.0)
DM, DC = 0.05, 5.0
DIFFUSION = {"u": DM, "q": DM, "v": DC, "w": DC, "p": DC}
M_BG = 14.0                 # total GEF  = integral of u + q + v + w over [0, L]
M_P = 3.444                 # total PAK  = integral of p + q over [0, L]; rho = M_P/M_BG = 0.246 (Table 2 rounds to 0.25)
N_GEF = M_BG / L            # the package takes mean concentrations: 2.0
N_PAK = M_P / L             # 0.492
TAU_C = 1.0                 # noise correlation time
N_GRID = 200
PERT = 0.01                 # relative amplitude of the initial perturbation


def trapz(y, x, axis=-1):
    """Trapezoid rule; the same arithmetic as numpy.trapz (deprecated in numpy 2)."""
    y = np.asarray(y)
    d = np.diff(x)
    shape = [1] * y.ndim
    shape[axis] = d.shape[0]
    d = d.reshape(shape)
    s1 = [slice(None)] * y.ndim
    s2 = [slice(None)] * y.ndim
    s1[axis] = slice(1, None)
    s2[axis] = slice(None, -1)
    return np.add.reduce(d * (y[tuple(s1)] + y[tuple(s2)]) / 2.0, axis)


def make_grid(n=N_GRID):
    return Grid(N=n, L=L * n / (n - 1), bc="neumann")


def trap_weights(n):
    dx = L / (n - 1)
    w = np.full(n, dx)
    w[0] = w[-1] = dx / 2
    return w


def hss():
    eqs = solve_homogeneous_equilibrium(N_GEF, N_PAK, PARAMS)
    if len(eqs) != 1:
        raise RuntimeError(f"expected one HSS, got {len(eqs)}: {eqs}")
    return eqs[0]


def initial_condition(rep, n=N_GRID):
    """HSS + 1 % relative Gaussian perturbation (seed = rep), perturbation mean removed per species."""
    eq = hss()
    rng = np.random.default_rng(rep)
    u0 = np.tile(eq[:, None], (1, n)) + PERT * eq[:, None] * rng.standard_normal((5, n))
    w = trap_weights(n)
    u0 = u0 - ((u0 * w).sum(1, keepdims=True) / L - eq[:, None])
    return u0


class MassFixSolver(MCRDSolver):
    """The package's Crank-Nicolson scheme (MCRDSolver._step) with one change: the trapezoid-weighted spatial
    mean of each noise increment (u and q separately) is subtracted before injection."""

    def run(self, u0, n_steps, *, save_every=1, noisy=False, noise_A=None,
            noise_tau_c=None, seed=None):
        u = u0.copy().astype(float)
        phi = reaction(u, self.params)
        N = self.grid.N
        w = trap_weights(N)
        rng = np.random.default_rng(seed)
        fields = {i: OUNoiseField(N=N, A=noise_A, tau_c=noise_tau_c, rng=rng) for i in (0, 1)}
        n_saved = n_steps // save_every + 1
        times = np.zeros(n_saved)
        snaps = np.zeros((n_saved, 5, N))
        snaps[0] = u
        k = 1
        for step in range(1, n_steps + 1):
            nr = np.zeros((5, N))
            for i, fld in fields.items():
                en = fld.value.copy()
                enp = fld.step(self.dt)
                inc = (self.dt / 2.0) * (en + enp)
                nr[i] = inc - (inc * w).sum() / L          # remove the spatial mean
            u, phi, it = self._step(u, phi, nr)
            if it == self.fp_maxiter:
                raise RuntimeError(f"fixed point did not converge at step {step}")
            if step % save_every == 0:
                times[k] = step * self.dt
                snaps[k] = u
                k += 1
        return SimulationResult(times=times, snapshots=snaps, mean=snaps.mean(2),
                                std=snaps.std(2), mass_gef=None, mass_pak=None,
                                grid=self.grid, noisy=True)


# ---------------------------------------------------------------- observables
def observables(times, snaps):
    """GEF = u + q + v + w, PAK = p + q; left half = [0, L/2].
    dG, dP: (left - right) half masses divided by M_BG.
    dn2_u = int_0^L cos(2 pi x / L) u dx (the shape axis of Figs. 3 and 4).
    dn2   = int_0^L cos(2 pi x / L) (u + q) dx / M_BG (a variant normalisation; not used in the manuscript)."""
    n = snaps.shape[2]
    x = np.linspace(0.0, L, n)
    g = snaps[:, 0] + snaps[:, 1] + snaps[:, 2] + snaps[:, 3]
    k = snaps[:, 4] + snaps[:, 1]
    uq = snaps[:, 0] + snaps[:, 1]
    xm = L / 2
    im = np.searchsorted(x, xm)            # first node right of the midpoint
    xl = np.r_[x[:im], xm]
    xr = np.r_[xm, x[im:]]

    def halves(y):
        ymid = y[:, im - 1] + (y[:, im] - y[:, im - 1]) * (xm - x[im - 1]) / (x[im] - x[im - 1])
        yl = np.concatenate([y[:, :im], ymid[:, None]], 1)
        yr = np.concatenate([ymid[:, None], y[:, im:]], 1)
        return trapz(yl, xl, axis=1), trapz(yr, xr, axis=1)

    gL, gR = halves(g)
    kL, kR = halves(k)
    dn2 = trapz(np.cos(2 * np.pi * x / L)[None, :] * uq, x, axis=1) / M_BG
    dn2_u = trapz(np.cos(2 * np.pi * x / L)[None, :] * snaps[:, 0], x, axis=1)
    u = snaps[:, 0]
    return dict(
        t=times,
        uL_max=u[:, :im].max(1), uR_max=u[:, im:].max(1),
        dG=(gL - gR) / M_BG, dP=(kL - kR) / M_BG, dn2=dn2, dn2_u=dn2_u,
        M_GEF=gL + gR, M_PAK=kL + kR,
        min_conc=snaps.min(axis=(1, 2)),
    )


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--out", required=True, help="output folder")
    ap.add_argument("--A", type=float, default=0.001, help="noise strength A = sigma^2 tau_c (0 = no noise)")
    ap.add_argument("--rep", type=int, default=0, help="replicate index (initial-condition and noise seeds)")
    ap.add_argument("--T", type=float, default=4000.0, help="end time (s)")
    ap.add_argument("--dt", type=float, default=0.0025, help="time step (s); 0.01 and 0.005 fail to converge")
    ap.add_argument("--save_dt", type=float, default=0.5, help="interval of the saved observables (s)")
    ap.add_argument("--snap_dt", type=float, default=2.0, help="interval of the saved u profiles (s)")
    ap.add_argument("--mode", choices=["massfix", "uncorrected"], default="massfix",
                    help="massfix: spatial mean of each noise increment removed, totals conserved (Fig. 4); uncorrected: noise added as generated, totals conserved only on average")
    ap.add_argument("--force", action="store_true", help="overwrite an existing output file")
    a = ap.parse_args(argv)

    os.makedirs(a.out, exist_ok=True)
    fn = os.path.join(a.out, f"A{a.A:g}_rep{a.rep}.npz")
    if os.path.exists(fn) and not a.force:
        print("exists, skip", fn)
        return
    grid = make_grid()
    u0 = initial_condition(a.rep)
    n_steps = int(round(a.T / a.dt))
    se = int(round(a.save_dt / a.dt))
    cls = MassFixSolver if (a.mode == "massfix" and a.A > 0) else MCRDSolver
    solver = cls(grid, DIFFUSION, PARAMS, dt=a.dt)
    seed = 100000 + int(round(a.A * 1e4)) * 100 + a.rep   # noise seed per (A, rep)
    t0 = time.time()
    status = "ok"
    try:
        if a.A > 0:
            res = solver.run(u0, n_steps, save_every=se, noisy=True, noise_A=a.A,
                             noise_tau_c=TAU_C, seed=seed)
        else:
            res = solver.run(u0, n_steps, save_every=se, noisy=False)
        times, snaps = res.times, res.snapshots
    except RuntimeError as exc:   # fixed-point non-convergence: recorded, not hidden
        status = f"failed: {exc}"
        times = snaps = None
    wall = time.time() - t0
    meta = dict(A=a.A, rep=a.rep, T=a.T, dt=a.dt, mode=a.mode, noise_seed=seed,
                ic_seed=a.rep, tau_c=TAU_C, N_grid=N_GRID, L=L, M_BG=M_BG, M_P=M_P,
                N_GEF=N_GEF, N_PAK=N_PAK, wall_s=wall, status=status)
    if times is None:
        np.savez(fn, meta=np.array(repr(meta)))
        print(fn, status, f"{wall:.1f}s")
        return
    obs = observables(times, snaps)
    ss = int(round(a.snap_dt / a.save_dt))
    np.savez_compressed(
        fn, meta=np.array(repr(meta)), **obs,
        snap_t=times[::ss], snap_u=snaps[::ss, 0].astype(np.float32),
        final_state=snaps[-1],
    )
    print(fn, status, f"{wall:.1f}s", "M_GEF drift %.3g M_PAK drift %.3g min %.3g" % (
        obs["M_GEF"][-1] - obs["M_GEF"][0], obs["M_PAK"][-1] - obs["M_PAK"][0], obs["min_conc"].min()))


if __name__ == "__main__":
    main()
