"""First-pass deliverable: does paper-style OU active noise change the dynamics of the
Gef/Pak1 MCRD system?

Runs one noiseless and one noisy trajectory from the same (equilibrium + small
perturbation) initial condition and the same random seed for the initial perturbation,
so the only difference between the two runs is the OU forcing itself. Produces:
  - a spacetime heatmap of each species for the noisy run
  - overlaid noisy-vs-noiseless spatial mean/std-dev traces per species
  - n_GEF(t) and n_pak(t) for both runs, to quantify how much the *unconstrained*
    noise on (u, q) drifts the two conserved masses

ALL PARAMETER VALUES BELOW ARE PLACEHOLDERS. They were chosen only to give a valid
positive homogeneous equilibrium and a numerically well-behaved run, so a first
execution is a mechanics check (does the machinery work, does it visibly show
something) -- not a physically calibrated result. Swap in real values before drawing
scientific conclusions.
"""

from __future__ import annotations

import numpy as np
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt

from gef_pak_noise.grid import Grid
from gef_pak_noise.model import ReactionParams, solve_homogeneous_equilibrium
from gef_pak_noise.solver import MCRDSolver, SimulationResult, _SPECIES

# ---------------------------------------------------------------------------
# CONFIG -- placeholders, see module docstring. Replace with fitted/known values.
# ---------------------------------------------------------------------------

REACTION_PARAMS = ReactionParams(a=1.0, b=1.0, c=1.0, d=1.0, e=1.0, f=1.0)  # PLACEHOLDER

# Membrane species (u, q) assumed slow; cytosolic (v, w, p) assumed fast -- see the
# diffusion-label discrepancy flagged in the plan. Flip this if it's backwards.
DIFFUSION = {"u": 0.1, "q": 0.1, "v": 5.0, "w": 5.0, "p": 5.0}  # PLACEHOLDER

N_GEF = 5.0  # PLACEHOLDER -- control parameter: total Gef
N_PAK = 3.0  # PLACEHOLDER -- control parameter: total Pak1
DOMAIN_LENGTH = 10.0  # PLACEHOLDER -- control parameter: domain length L
BOUNDARY_CONDITION = "periodic"  # flagged assumption -- switch to "neumann" if needed

N_GRID = 200
DT = 0.01
T_FINAL = 20.0
SAVE_EVERY = 20

NOISE_A = 0.05  # PLACEHOLDER -- integrated noise strength (paper's Eq. 3)
NOISE_TAU_C = 1.0  # PLACEHOLDER -- OU correlation time
PERTURBATION_AMPLITUDE = 0.01  # relative amplitude of initial random perturbation
SEED = 0

OUTPUT_DIR = "outputs"


def build_initial_condition(equilibrium: np.ndarray, N: int, rng: np.random.Generator) -> np.ndarray:
    """Equilibrium + small relative random perturbation, same for noisy/noiseless runs."""
    u0 = np.tile(equilibrium[:, None], (1, N)).astype(float)
    perturbation = PERTURBATION_AMPLITUDE * equilibrium[:, None] * rng.standard_normal((5, N))
    return u0 + perturbation


def main():
    import os

    os.makedirs(OUTPUT_DIR, exist_ok=True)

    equilibria = solve_homogeneous_equilibrium(N_GEF, N_PAK, REACTION_PARAMS)
    if not equilibria:
        raise RuntimeError(
            "No valid positive homogeneous equilibrium found for the current "
            "REACTION_PARAMS / N_GEF / N_PAK -- these placeholders need adjusting "
            "before a run is possible."
        )
    print(f"Found {len(equilibria)} valid equilibrium/equilibria:")
    for eq in equilibria:
        print("  ", dict(zip(_SPECIES, eq)))
    equilibrium = equilibria[0]
    print(f"Using equilibrium: {dict(zip(_SPECIES, equilibrium))}")

    grid = Grid(N=N_GRID, L=DOMAIN_LENGTH, bc=BOUNDARY_CONDITION)
    print(f"Boundary condition: {BOUNDARY_CONDITION} (flagged assumption, see plan)")

    n_steps = int(round(T_FINAL / DT))

    rng_ic = np.random.default_rng(SEED)
    u0 = build_initial_condition(equilibrium, grid.N, rng_ic)

    solver = MCRDSolver(grid, DIFFUSION, REACTION_PARAMS, dt=DT)

    print("Running noiseless trajectory...")
    result_clean = solver.run(u0, n_steps, save_every=SAVE_EVERY, noisy=False)

    print("Running noisy trajectory...")
    result_noisy = solver.run(
        u0,
        n_steps,
        save_every=SAVE_EVERY,
        noisy=True,
        noise_A=NOISE_A,
        noise_tau_c=NOISE_TAU_C,
        seed=SEED + 1,  # noise realization seed; initial condition (u0) is identical
    )

    make_plots(result_clean, result_noisy, equilibrium)
    save_arrays(result_clean, result_noisy)

    # Quantify the two things this first pass is meant to test.
    report(result_clean, result_noisy, equilibrium)


def make_plots(clean: SimulationResult, noisy: SimulationResult, equilibrium: np.ndarray):
    x = clean.grid.x

    # 1. Spacetime heatmaps for the noisy run, one panel per species.
    fig, axes = plt.subplots(1, 5, figsize=(20, 4), sharey=True)
    for i, ax in enumerate(axes):
        im = ax.pcolormesh(
            x, noisy.times, noisy.snapshots[:, i, :], shading="auto", cmap="magma"
        )
        ax.set_title(_SPECIES[i])
        ax.set_xlabel("position")
        fig.colorbar(im, ax=ax, fraction=0.046)
    axes[0].set_ylabel("time")
    fig.suptitle("Noisy run: spacetime concentration profiles")
    fig.tight_layout()
    fig.savefig(f"{OUTPUT_DIR}/spacetime_noisy.png", dpi=150)
    plt.close(fig)

    # 2. Noisy-vs-noiseless mean and std-dev per species.
    fig, axes = plt.subplots(2, 5, figsize=(20, 6), sharex=True)
    for i in range(5):
        axes[0, i].plot(clean.times, clean.mean[:, i], label="noiseless")
        axes[0, i].plot(noisy.times, noisy.mean[:, i], label="noisy", alpha=0.8)
        axes[0, i].axhline(equilibrium[i], color="gray", linestyle=":", label="equilibrium")
        axes[0, i].set_title(f"{_SPECIES[i]}: spatial mean")

        axes[1, i].plot(clean.times, clean.std[:, i], label="noiseless")
        axes[1, i].plot(noisy.times, noisy.std[:, i], label="noisy", alpha=0.8)
        axes[1, i].set_title(f"{_SPECIES[i]}: spatial std-dev")
        axes[1, i].set_xlabel("time")
    axes[0, 0].legend(fontsize=8)
    fig.tight_layout()
    fig.savefig(f"{OUTPUT_DIR}/mean_std_comparison.png", dpi=150)
    plt.close(fig)

    # 3. Conserved-mass drift.
    fig, axes = plt.subplots(1, 2, figsize=(10, 4))
    axes[0].plot(clean.times, clean.mass_gef, label="noiseless")
    axes[0].plot(noisy.times, noisy.mass_gef, label="noisy")
    axes[0].axhline(N_GEF, color="gray", linestyle=":", label="target N_GEF")
    axes[0].set_title("n_GEF(t) = <u+q+v+w>")
    axes[0].set_xlabel("time")
    axes[0].legend(fontsize=8)

    axes[1].plot(clean.times, clean.mass_pak, label="noiseless")
    axes[1].plot(noisy.times, noisy.mass_pak, label="noisy")
    axes[1].axhline(N_PAK, color="gray", linestyle=":", label="target N_pak")
    axes[1].set_title("n_pak(t) = <p+q>")
    axes[1].set_xlabel("time")
    axes[1].legend(fontsize=8)
    fig.tight_layout()
    fig.savefig(f"{OUTPUT_DIR}/mass_drift.png", dpi=150)
    plt.close(fig)


def save_arrays(clean: SimulationResult, noisy: SimulationResult):
    np.savez(
        f"{OUTPUT_DIR}/trajectories.npz",
        times=clean.times,
        clean_snapshots=clean.snapshots,
        noisy_snapshots=noisy.snapshots,
        clean_mean=clean.mean,
        noisy_mean=noisy.mean,
        clean_std=clean.std,
        noisy_std=noisy.std,
        clean_mass_gef=clean.mass_gef,
        noisy_mass_gef=noisy.mass_gef,
        clean_mass_pak=clean.mass_pak,
        noisy_mass_pak=noisy.mass_pak,
        species=_SPECIES,
    )


def report(clean: SimulationResult, noisy: SimulationResult, equilibrium: np.ndarray):
    print("\n--- Summary ---")
    # Compare against each run's own initial mass, not the nominal N_GEF/N_pak target:
    # the random initial perturbation itself has a finite-sample mean that generally
    # isn't exactly N_GEF/N_pak, and the solver (correctly) conserves whatever mass the
    # initial condition actually has. Comparing to the nominal target would conflate
    # that harmless sampling offset with actual solver/noise-driven drift.
    print(f"Noiseless run: mass_gef drifted {clean.mass_gef[-1] - clean.mass_gef[0]:+.3e} "
          f"from its initial value ({clean.mass_gef[0]:.6f}) -- should be ~0 to solver "
          f"tolerance, since reaction+diffusion conserve mass exactly by construction.")
    print(f"Noisy run:     mass_gef drifted {noisy.mass_gef[-1] - noisy.mass_gef[0]:+.3e} "
          f"from its initial value ({noisy.mass_gef[0]:.6f}) -- nonzero is expected, "
          f"since noise on (u, q) here is deliberately unconstrained.")
    print(f"Noiseless run: mass_pak drifted {clean.mass_pak[-1] - clean.mass_pak[0]:+.3e}")
    print(f"Noisy run:     mass_pak drifted {noisy.mass_pak[-1] - noisy.mass_pak[0]:+.3e}")

    print("\nFinal spatial std-dev per species "
          f"({', '.join(_SPECIES)}):")
    print(f"  noiseless: {np.array2string(clean.std[-1], precision=4)}")
    print(f"  noisy:     {np.array2string(noisy.std[-1], precision=4)}")
    print("(the noiseless run's std-dev decaying toward ~0 means the initial "
          "perturbation simply relaxes back to the homogeneous equilibrium without "
          "noise -- i.e. this parameter set has no pattern-forming instability of its "
          "own; the noisy run's std-dev settling at a nonzero plateau shows the OU "
          "forcing sustains persistent spatial heterogeneity that wouldn't otherwise "
          "exist. A ratio between the two isn't reported since it blows up "
          "uninformatively once the noiseless value is near machine precision.)")


if __name__ == "__main__":
    main()
