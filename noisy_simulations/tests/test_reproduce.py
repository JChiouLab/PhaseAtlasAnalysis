"""Reproduction test: rerun the first 40 s of two Fig. 4 runs and compare with the stored production output.

Runs: A = 0.001, rep 0 (noise, mode massfix) and A = 0, rep 0 (the noise-free control). The reference file
reference_first40s.npz holds the first 40 s of the original output files, cut from them (not from a rerun).
On the original environment (Python 3.11.8, numpy 1.26.4, scipy 1.12.0) the agreement is exact.
Other versions may differ at round-off; the test passes below a tolerance of 1e-9 and prints the largest difference.

Run: python3 tests/test_reproduce.py   (about 10 s)
"""
import os
import subprocess
import sys
import tempfile

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
DRIVER = os.path.join(HERE, "..", "noise_driver.py")
REF = os.path.join(HERE, "reference_first40s.npz")
KEYS = ["t", "uL_max", "uR_max", "dG", "dP", "dn2", "dn2_u", "M_GEF", "M_PAK", "min_conc", "snap_t", "snap_u"]
TOL = 1e-9


def main():
    ref = np.load(REF)
    worst = 0.0
    with tempfile.TemporaryDirectory() as out:
        for A, tag in [(0.001, "A0.001_rep0"), (0.0, "A0_rep0")]:
            subprocess.run([sys.executable, DRIVER, "--out", out, "--A", str(A), "--rep", "0", "--T", "40"],
                           check=True)
            z = np.load(os.path.join(out, f"{tag}.npz"))
            for k in KEYS:
                r = ref[f"{tag}__{k}"]
                d = float(np.max(np.abs(z[k].astype(float) - r.astype(float))))
                worst = max(worst, d)
                print(f"{tag:12s} {k:9s} max |diff| = {d:.3g}")
    print(f"largest difference {worst:.3g}")
    if worst > TOL:
        print("FAIL")
        sys.exit(1)
    print("PASS" + (" (bit-exact)" if worst == 0 else ""))


if __name__ == "__main__":
    main()
