#!/usr/bin/env python3
"""
emissivity_rk45_tol_sweep.py

Sweeps the RK45 DOPRI5 error tolerance (rk45_tol) over a range of values and
plots how the photon-sphere separatrix discrepancy and wall-clock runtime vary.

Calls:  bin/emissivity_rk45_plot <csv_path> <tol>
Reads the CSV, computes per-bin emissivity deviation  (RK45 - RK4) / RK4,
and records:
  - RMS of the emissivity deviation across all populated bins
  - Max absolute emissivity deviation
  - RK45 wall-clock time (sum of t column over all rays, as a proxy; the
    binary doesn't report wall time for the sweep, so we time the subprocess)

Output: emissivity_rk45_tol_sweep.png  (3-panel figure)
"""

import subprocess
import tempfile
import os
import time
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

# ---------------------------------------------------------------------------
# Configuration
# ---------------------------------------------------------------------------
PROJECT_ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), "../.."))
BIN          = os.path.join(PROJECT_ROOT, "bin", "emissivity_rk45_plot")
CMAKE        = "/opt/local/bin/cmake"
BUILD_DIR    = os.path.join(PROJECT_ROOT, "build")
OUTPUT_PNG   = os.path.join(PROJECT_ROOT, "emissivity_rk45_tol_sweep.png")

# Tolerances to sweep (coarse → tight)
TOL_VALUES = [1e-6, 3e-7, 1e-7, 3e-8, 1e-8, 3e-9, 1e-9, 3e-10, 1e-10]

# Minimum ray count in a bin for it to be included in the statistics
MIN_BIN_RAYS = 50

# ---------------------------------------------------------------------------
def build_binary():
    print("Building emissivity_rk45_plot...")
    result = subprocess.run(
        [CMAKE, "--build", BUILD_DIR, "--target", "emissivity_rk45_plot"],
        capture_output=True, text=True
    )
    if result.returncode != 0:
        print("Build FAILED:\n", result.stderr)
        raise RuntimeError("Build failed")
    print("Build OK.\n")


def run_one(tol):
    """Run the binary for the given tolerance, return (df, wall_seconds)."""
    with tempfile.NamedTemporaryFile(suffix=".csv", delete=False) as f:
        csv_path = f.name
    try:
        t0 = time.perf_counter()
        result = subprocess.run(
            [BIN, csv_path, str(tol)],
            capture_output=True, text=True
        )
        wall = time.perf_counter() - t0
        if result.returncode != 0:
            raise RuntimeError(f"Binary failed for tol={tol}:\n{result.stderr}")

        # Read CSV, skipping comment lines (# rk45_tol=...)
        rows = []
        header = None
        with open(csv_path) as fh:
            for line in fh:
                line = line.strip()
                if not line or line.startswith("#"):
                    continue
                if header is None:
                    header = line.split()
                    continue
                rows.append([float(x) for x in line.split()])

        if not rows:
            return None, wall

        data = {col: np.array([r[i] for r in rows])
                for i, col in enumerate(header)}
        return data, wall
    finally:
        try:
            os.unlink(csv_path)
        except OSError:
            pass


# ---------------------------------------------------------------------------
def compute_metrics(data):
    """Return (rms_dev, max_dev) of emissivity deviation for populated bins."""
    mask = (data["N_rk4"] >= MIN_BIN_RAYS) & (data["N_rk45"] >= MIN_BIN_RAYS)
    e4  = data["emis_rk4"][mask]
    e45 = data["emis_rk45"][mask]
    # guard against zero
    valid = e4 > 0
    rel   = (e45[valid] - e4[valid]) / e4[valid]
    if len(rel) == 0:
        return np.nan, np.nan
    return float(np.sqrt(np.mean(rel**2))), float(np.max(np.abs(rel)))


# ---------------------------------------------------------------------------
def main():
    build_binary()

    # Run RK4 once to get baseline (tol argument ignored for RK4 run)
    print(f"Running for tol = 1e-8 (baseline)...")
    baseline_data, _ = run_one(1e-8)

    tols       = []
    rms_devs   = []
    max_devs   = []
    wall_times = []

    for tol in TOL_VALUES:
        print(f"  tol = {tol:.1e} ... ", end="", flush=True)
        data, wall = run_one(tol)
        if data is None:
            print("no data, skipping")
            continue
        rms, mx = compute_metrics(data)
        tols.append(tol)
        rms_devs.append(rms)
        max_devs.append(mx)
        wall_times.append(wall)
        print(f"done  RMS dev = {rms:.3f}  max dev = {mx:.3f}  wall = {wall:.1f}s")

    tols       = np.array(tols)
    rms_devs   = np.array(rms_devs)
    max_devs   = np.array(max_devs)
    wall_times = np.array(wall_times)

    # Default tolerance index for reference line
    default_tol = 1e-8
    ref_idx = np.argmin(np.abs(tols - default_tol)) if len(tols) else None

    # -----------------------------------------------------------------------
    fig, axes = plt.subplots(3, 1, figsize=(8, 10), sharex=True)
    fig.suptitle(
        "RK45 (DOPRI5) tolerance sweep\n"
        r"spin = 0.998, lamppost $r=5\,r_g$",
        fontsize=13
    )

    # Panel 1 — RMS emissivity deviation
    axes[0].plot(tols, rms_devs * 100, "o-", color="C0")
    axes[0].set_ylabel("RMS emissivity deviation (%)")
    axes[0].set_xscale("log")
    axes[0].invert_xaxis()
    axes[0].grid(True, which="both", alpha=0.3)
    if ref_idx is not None:
        axes[0].axvline(default_tol, color="k", ls="--", lw=0.8, label=f"default tol={default_tol:.0e}")
        axes[0].legend(fontsize=9)

    # Panel 2 — max emissivity deviation
    axes[1].plot(tols, max_devs * 100, "s-", color="C1")
    axes[1].set_ylabel("Max |emissivity deviation| (%)")
    axes[1].set_xscale("log")
    axes[1].grid(True, which="both", alpha=0.3)
    if ref_idx is not None:
        axes[1].axvline(default_tol, color="k", ls="--", lw=0.8)

    # Panel 3 — wall-clock time
    axes[2].plot(tols, wall_times, "^-", color="C2")
    axes[2].set_ylabel("Wall-clock time (s)\n(RK4 + RK45)")
    axes[2].set_xlabel("RK45 error tolerance (tighter →)")
    axes[2].set_xscale("log")
    axes[2].grid(True, which="both", alpha=0.3)
    if ref_idx is not None:
        axes[2].axvline(default_tol, color="k", ls="--", lw=0.8)

    plt.tight_layout()
    plt.savefig(OUTPUT_PNG, dpi=150)
    print(f"\nSaved: {OUTPUT_PNG}")
    plt.show()


if __name__ == "__main__":
    main()
