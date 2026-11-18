"""Ornstein-Uhlenbeck active noise field, following the paper's Eqs. (2)-(4) and (A1).

The paper parameterizes the noise by the integrated noise strength
    A = sigma_eta**2 * tau_c
(Eq. 3) and treats A, rather than sigma_eta, as the primary control parameter, since
A is what they show actually controls the transition between dynamical regimes.

The noise field's own SDE (Eq. 4):
    d eta / dt = -eta / tau_c + sqrt(2 A / tau_c) * xi(l, t)
is discretized with Euler-Maruyama (Eq. A1):
    eta_i(t + dt) = (1 - dt/tau_c) * eta_i(t) + sqrt(2*A*dt) / tau_c * xi_i(t)
where xi_i(t) are iid standard Gaussians, independent across grid points i and
across the two noisy species (u, q) here.

Numerical implementation note (flagging an interpretation rather than resolving it
silently): the paper's Appendix A is terse about exactly when the noise increment is
applied relative to the Crank-Nicolson sub-step for the deterministic reaction-diffusion
terms. Here, eta is advanced once per step *before* the fixed-point iteration (since it
does not depend on u/v/... at all), and enters the Crank-Nicolson right-hand side with
trapezoidal averaging (eta^n + eta^{n+1})/2, consistent with how the reaction term phi
is treated in Eqs. (A5)-(A6). If exact fidelity to the reference implementation's
operator-splitting order matters, this is the piece to revisit.
"""

from dataclasses import dataclass

import numpy as np


@dataclass
class OUNoiseField:
    """One independent Ornstein-Uhlenbeck noise field on an N-point spatial grid."""

    N: int
    A: float
    tau_c: float
    rng: np.random.Generator
    value: np.ndarray = None  # type: ignore[assignment]

    def __post_init__(self):
        if self.value is None:
            # Initialize from the process's own stationary distribution,
            # Var(eta) = A / tau_c = sigma_eta**2, rather than from eta=0, so the
            # noise doesn't have an artificial "spin-up" transient at t=0.
            sigma_eta = np.sqrt(self.A / self.tau_c)
            self.value = sigma_eta * self.rng.standard_normal(self.N)

    def step(self, dt: float) -> np.ndarray:
        """Advance eta by one Euler-Maruyama step (Eq. A1) and return the new value."""
        if self.tau_c <= 0:
            # White-noise limit: A = sigma_eta**2 * tau_c with tau_c -> 0 at fixed A
            # diverges (Eq. 4 is singular); if you need this limit, sample white
            # noise directly instead of routing it through this OU process.
            raise ValueError(
                "tau_c must be > 0; the tau_c -> 0 white-noise limit is singular "
                "for the OU parameterization (see Eq. 4) and needs separate handling."
            )
        xi = self.rng.standard_normal(self.N)
        self.value = (1.0 - dt / self.tau_c) * self.value + (
            np.sqrt(2.0 * self.A * dt) / self.tau_c
        ) * xi
        return self.value
