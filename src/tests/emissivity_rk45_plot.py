#!/usr/bin/env python3
"""
emissivity_rk45_plot.py

Builds and runs emissivity_rk45_plot (C++ binary), then plots the per-bin
ray count, emissivity, mean redshift, and mean arrival time for both the
RK4 and RK45 integrators.

Run from the project root:
    python3 src/tests/emissivity_rk45_plot.py
"""

import subprocess
import tempfile
import os
import sys

import numpy as np
import matplotlib.pyplot as plt

# ---------------------------------------------------------------------------
# Locate project root and binary
# ---------------------------------------------------------------------------
SCRIPT_DIR  = os.path.dirname(os.path.abspath(__file__))
PROJECT_ROOT = os.path.abspath(os.path.join(SCRIPT_DIR, "..", ".."))
BUILD_DIR    = os.path.join(PROJECT_ROOT, "build")
BINARY       = os.path.join(PROJECT_ROOT, "bin", "emissivity_rk45_plot")
CMAKE        = "/opt/local/bin/cmake"

# ---------------------------------------------------------------------------
# Build the binary
# ---------------------------------------------------------------------------
print("Building emissivity_rk45_plot...")
result = subprocess.run(
    [CMAKE, "--build", BUILD_DIR, "--target", "emissivity_rk45_plot"],
    capture_output=True, text=True
)
if result.returncode != 0:
    print("Build failed:")
    print(result.stdout)
    print(result.stderr)
    sys.exit(1)
print("Build OK.")

# ---------------------------------------------------------------------------
# Run the binary, writing CSV to a temp file
# ---------------------------------------------------------------------------
with tempfile.NamedTemporaryFile(suffix=".csv", delete=False) as tmp:
    csv_path = tmp.name

try:
    print("Running integrators (this takes ~10 s)...")
    result = subprocess.run([BINARY, csv_path], capture_output=False, text=True)
    if result.returncode != 0:
        print("Binary failed with exit code", result.returncode)
        sys.exit(1)

    # -----------------------------------------------------------------------
    # Load CSV
    # -----------------------------------------------------------------------
    data = np.loadtxt(csv_path, skiprows=1)   # skip header line
finally:
    os.unlink(csv_path)

r         = data[:, 0]
N_rk4     = data[:, 1]
N_rk45    = data[:, 2]
emis_rk4  = data[:, 3]
emis_rk45 = data[:, 4]
redsh_rk4 = data[:, 5]
redsh_rk45= data[:, 6]
time_rk4  = data[:, 7]
time_rk45 = data[:, 8]

# Zero emissivity/redshift/time values mean no rays in that bin — mask them
def masked(arr):
    a = arr.copy().astype(float)
    a[a == 0] = np.nan
    return a

# ---------------------------------------------------------------------------
# Plot
# ---------------------------------------------------------------------------
fig, axes = plt.subplots(2, 2, figsize=(11, 8))
fig.suptitle(
    "RK4 vs RK45 emissivity comparison\n"
    r"spin = 0.998,  lamppost $r=5\,r_g$,  $\theta=10^{-3}$ rad",
    fontsize=12
)

kw4  = dict(color="steelblue",  lw=1.8, marker="o", ms=4, label="RK4")
kw45 = dict(color="tomato",     lw=1.8, marker="s", ms=4, label="RK45",
            linestyle="--")

# --- Ray count ---
ax = axes[0, 0]
ax.plot(r, N_rk4,  **kw4)
ax.plot(r, N_rk45, **kw45)
ax.set_xscale("log")
ax.set_xlabel(r"$r\ [r_g]$")
ax.set_ylabel("Rays per bin")
ax.set_title("Ray count")
ax.legend()
ax.grid(True, which="both", ls=":", alpha=0.5)

# --- Emissivity ---
ax = axes[0, 1]
ax.plot(r, masked(emis_rk4),  **kw4)
ax.plot(r, masked(emis_rk45), **kw45)
ax.set_xscale("log")
ax.set_yscale("log")
ax.set_xlabel(r"$r\ [r_g]$")
ax.set_ylabel(r"Emissivity $[\mathrm{arb.}]$")
ax.set_title("Disc emissivity")
ax.legend()
ax.grid(True, which="both", ls=":", alpha=0.5)

# --- Mean redshift ---
ax = axes[1, 0]
ax.plot(r, masked(redsh_rk4),  **kw4)
ax.plot(r, masked(redsh_rk45), **kw45)
ax.set_xscale("log")
ax.set_xlabel(r"$r\ [r_g]$")
ax.set_ylabel(r"Mean redshift $g$")
ax.set_title("Mean photon redshift")
ax.legend()
ax.grid(True, which="both", ls=":", alpha=0.5)

# --- Mean arrival time ---
ax = axes[1, 1]
ax.plot(r, masked(time_rk4),  **kw4)
ax.plot(r, masked(time_rk45), **kw45)
ax.set_xscale("log")
ax.set_xlabel(r"$r\ [r_g]$")
ax.set_ylabel(r"Mean arrival time $[r_g/c]$")
ax.set_title("Mean photon arrival time")
ax.legend()
ax.grid(True, which="both", ls=":", alpha=0.5)

plt.tight_layout()

outfile = os.path.join(PROJECT_ROOT, "emissivity_rk45_comparison.png")
plt.savefig(outfile, dpi=150)
print(f"Plot saved to {outfile}")
plt.show()
