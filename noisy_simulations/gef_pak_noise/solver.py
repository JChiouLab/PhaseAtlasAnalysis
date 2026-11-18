"""Crank-Nicolson time-stepper for the 5-species Gef/Pak1 MCRD system, generalizing
Appendix A of the paper (Eqs. A5-A14) from 2 species to 5.

Because diffusion is diagonal (D = Diag(D_u, D_q, D_v, D_w, D_p), no cross-diffusion),
the linear (diffusive) part of the Crank-Nicolson update decouples into 5 independent
per-species linear solves -- exactly the structure of the paper's Au/Av, Bu/Bv blocks,
just with 5 blocks instead of 2. The nonlinear reaction term f(u_vec) couples the 5
species pointwise at each grid location (not across grid locations), so it's handled
with the same fixed-point (Picard) iteration the paper uses for its nonlinear phi term,
generalized to vector-valued f.

Noise (independent OU fields on u and q, from noise.py) is injected as an explicit
source term in the Crank-Nicolson right-hand side, trapezoidal-averaged like every
other explicit term in the scheme -- see the note at the top of noise.py.
"""

from dataclasses import dataclass, field

import numpy as np
import scipy.sparse as sp
import scipy.sparse.linalg as spla

from .grid import Grid, build_laplacian
from .model import ReactionParams, reaction
from .noise import OUNoiseField

_SPECIES = ("u", "q", "v", "w", "p")
_NOISY_SPECIES_IDX = (0, 1)  # u, q


@dataclass
class SimulationResult:
    times: np.ndarray
    snapshots: np.ndarray  # shape (n_saved, 5, N)
    mean: np.ndarray  # shape (n_saved, 5) -- spatial mean per species over time
    std: np.ndarray  # shape (n_saved, 5) -- spatial std-dev per species over time
    mass_gef: np.ndarray  # shape (n_saved,) -- spatial average of u+q+v+w
    mass_pak: np.ndarray  # shape (n_saved,) -- spatial average of p+q
    grid: Grid
    noisy: bool


class MCRDSolver:
    def __init__(
        self,
        grid: Grid,
        D: dict,
        params: ReactionParams,
        dt: float,
        fp_tol: float = 1e-9,
        fp_maxiter: int = 50,
    ):
        """
        Parameters
        ----------
        grid : Grid
        D : dict mapping species name -> diffusion constant, e.g.
            {"u": D_u, "q": D_q, "v": D_v, "w": D_w, "p": D_p}
        params : ReactionParams
        dt : float, timestep
        """
        self.grid = grid
        self.params = params
        self.dt = dt
        self.fp_tol = fp_tol
        self.fp_maxiter = fp_maxiter

        missing = [s for s in _SPECIES if s not in D]
        if missing:
            raise ValueError(f"Missing diffusion constants for species: {missing}")
        self.D = np.array([D[s] for s in _SPECIES])

        L = build_laplacian(grid)
        I = sp.identity(grid.N, format="csr")

        self._factorized = []
        self._B = []
        for D_s in self.D:
            lam = D_s * dt / 2.0
            A_s = (I - lam * L).tocsc()
            B_s = I + lam * L
            self._factorized.append(spla.splu(A_s))
            self._B.append(B_s)

    def _step(
        self, u_vec: np.ndarray, phi_n: np.ndarray, noise_rhs: np.ndarray | None
    ) -> tuple[np.ndarray, np.ndarray, int]:
        """One Crank-Nicolson step with fixed-point iteration on the reaction term.

        Returns (u_next, phi_next, n_iterations_used).
        """
        dt = self.dt
        base_rhs = np.stack(
            [self._B[s] @ u_vec[s] for s in range(5)], axis=0
        )
        if noise_rhs is not None:
            base_rhs = base_rhs + noise_rhs

        phi_k = phi_n.copy()
        u_new = u_vec.copy()
        for it in range(1, self.fp_maxiter + 1):
            rhs = base_rhs + (dt / 2.0) * (phi_n + phi_k)
            for s in range(5):
                u_new[s] = self._factorized[s].solve(rhs[s])
            phi_new = reaction(u_new, self.params)
            resid = np.max(np.abs(phi_new - phi_k))
            phi_k = phi_new
            if resid < self.fp_tol:
                return u_new, phi_new, it

        return u_new, phi_k, self.fp_maxiter

    def run(
        self,
        u0: np.ndarray,
        n_steps: int,
        *,
        save_every: int = 1,
        noisy: bool = False,
        noise_A: float | None = None,
        noise_tau_c: float | None = None,
        seed: int | None = None,
    ) -> SimulationResult:
        """Integrate n_steps forward from initial condition u0 (shape (5, N))."""
        u = u0.copy().astype(float)
        phi = reaction(u, self.params)
        N = self.grid.N

        rng = np.random.default_rng(seed)
        noise_fields = None
        if noisy:
            if noise_A is None or noise_tau_c is None:
                raise ValueError("noise_A and noise_tau_c are required when noisy=True")
            noise_fields = {
                idx: OUNoiseField(N=N, A=noise_A, tau_c=noise_tau_c, rng=rng)
                for idx in _NOISY_SPECIES_IDX
            }

        n_saved = n_steps // save_every + 1
        times = np.zeros(n_saved)
        snapshots = np.zeros((n_saved, 5, N))
        mean = np.zeros((n_saved, 5))
        std = np.zeros((n_saved, 5))
        mass_gef = np.zeros(n_saved)
        mass_pak = np.zeros(n_saved)

        def record(idx, t, state):
            times[idx] = t
            snapshots[idx] = state
            mean[idx] = state.mean(axis=1)
            std[idx] = state.std(axis=1)
            mass_gef[idx] = (state[0] + state[1] + state[2] + state[3]).mean()
            mass_pak[idx] = (state[4] + state[1]).mean()

        record(0, 0.0, u)
        save_idx = 1

        for step in range(1, n_steps + 1):
            noise_rhs = None
            if noisy:
                noise_rhs = np.zeros((5, N))
                for idx, field_ in noise_fields.items():
                    eta_n = field_.value.copy()
                    eta_np1 = field_.step(self.dt)
                    noise_rhs[idx] = (self.dt / 2.0) * (eta_n + eta_np1)

            u, phi, n_iter = self._step(u, phi, noise_rhs)

            if n_iter == self.fp_maxiter:
                raise RuntimeError(
                    f"Fixed-point iteration did not converge by step {step} "
                    f"(residual tolerance {self.fp_tol}); reduce dt or check parameters."
                )

            if step % save_every == 0:
                record(save_idx, step * self.dt, u)
                save_idx += 1

        return SimulationResult(
            times=times,
            snapshots=snapshots,
            mean=mean,
            std=std,
            mass_gef=mass_gef,
            mass_pak=mass_pak,
            grid=self.grid,
            noisy=noisy,
        )
