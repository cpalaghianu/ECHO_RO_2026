#!/usr/bin/env python3
"""Sensitivity of the ECHO-LUCAS comparison to the position (origin) of the grid (modifiable areal unit problem, zoning effect).
For each cell size (20, 30, 50 km) the Stereo 70 grid is shifted over 5 x 5 offsets spanning one cell; main analysis, 100 x 100 km blocks,
bootstrap reduced to 100 replicates (only point estimates are summarised)."""
import subprocess, json, os, itertools, sys
import pandas as pd
HERE = os.path.dirname(os.path.abspath(__file__)); ROOT = os.path.join(HERE, "..")
# ECHO export (pH is the only variable used) and the LUCAS 2018 topsoil CSV for Romania (ESDAC, access-restricted, not in this repository)
E_ = sys.argv[1] if len(sys.argv) > 1 else os.path.join(ROOT, "data", "input", "samples_RO_FINAL.csv")
L_ = sys.argv[2] if len(sys.argv) > 2 else os.path.join(ROOT, "data", "external", "LUCAS-SOIL-2018 RO.csv")
BASE = (134140.3191, 753166.4394)
rows = []
out = "sens_origin"; os.makedirs(out, exist_ok=True)
for cell in (20, 30, 50):
    step = cell / 5
    for i, j in itertools.product(range(5), range(5)):
        ox, oy = BASE[0] + i * step * 1000, BASE[1] + j * step * 1000
        d = f"{out}/c{cell}_{i}{j}"
        subprocess.run([sys.executable, os.path.join(HERE, "lucas_compare.py"), "--echo", E_, "--lucas", L_, "--out", d, "--cell", str(cell * 1000), "--crs", "3844",
                        "--origin", f"{ox},{oy}", "--boot", "100"], stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL, check=True)
        a = json.load(open(d + "/lucas_summary.json")); m = a["main"]
        rows.append({"cell_km": cell, "i": i, "j": j, "n_cells": m["n_cells"], "rho": m["rho"]["est"], "kappa": m["kappa"]["est"],
                     "same_pct": m["same"]["est"], "meddiff": m["meddiff"]["est"], "rho_within100": m["rho_within"]["est"]})
df = pd.DataFrame(rows); df.to_csv(f"{out}/origin_robustness.csv", index=False)
s = df.groupby("cell_km").agg(n_cells_min=("n_cells", "min"), n_cells_max=("n_cells", "max"),
                              rho_min=("rho", "min"), rho_med=("rho", "median"), rho_max=("rho", "max"),
                              kappa_min=("kappa", "min"), kappa_med=("kappa", "median"), kappa_max=("kappa", "max"),
                              within_min=("rho_within100", "min"), within_med=("rho_within100", "median"), within_max=("rho_within100", "max"))
s.to_csv(f"{out}/origin_robustness_summary.csv"); print(s.round(2).to_string())
