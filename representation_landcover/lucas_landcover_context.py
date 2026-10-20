#!/usr/bin/env python3
"""Context check for the ECHO land-cover contrasts: the same contrasts in the independent LUCAS 2018 Topsoil survey (Romania, 0-20 cm,
laboratory pH in water and organic carbon). LUCAS land cover = first letter of LC (B cropland, C woodland, E grassland).
Descriptive only: different points, depth, year and land-cover classification than ECHO/CORINE."""
import sys, numpy as np, pandas as pd
from scipy import stats
src = sys.argv[1]; out = sys.argv[2]
d = pd.read_csv(src, low_memory=False)
if "NUTS_0" in d: d = d[d.NUTS_0.astype(str) == "RO"]      # one record in the file is labelled as Italy
for c in ["pH_H2O", "OC"]:
    d[c] = pd.to_numeric(d[c].astype(str).str.replace(",", ".", regex=False), errors="coerce")
d["lc"] = d.LC.astype(str).str[0]
rng = np.random.default_rng(20261019)
def delta(a, b):
    a, b = a.dropna().values, b.dropna().values
    return 2 * stats.mannwhitneyu(a, b).statistic / (len(a) * len(b)) - 1
rows = []
for name, g1, g2 in [("Woodland (C) vs cropland (B)", "C", "B"), ("Grassland (E) vs cropland (B)", "E", "B")]:
    for v, lab in [("OC", "Organic carbon (g/kg)"), ("pH_H2O", "pH (H2O)")]:
        a, b = d.loc[d.lc == g1, v].dropna(), d.loc[d.lc == g2, v].dropna()
        bs = [delta(pd.Series(rng.choice(a.values, len(a))), pd.Series(rng.choice(b.values, len(b)))) for _ in range(4000)]
        rows.append({"contrast": name, "variable": lab, "n1": len(a), "n2": len(b), "median1": a.median(), "median2": b.median(),
                     "cliffs_delta": delta(a, b), "ci_lo": np.percentile(bs, 2.5), "ci_hi": np.percentile(bs, 97.5),
                     "p": stats.mannwhitneyu(a, b).pvalue})
r = pd.DataFrame(rows); r.to_csv(out + "/LUCAS_RO_landcover_contrasts.csv", index=False)
m = d[["OC", "pH_H2O"]].dropna()
rho = stats.spearmanr(m.OC, m.pH_H2O)
open(out + "/LUCAS_RO_OC_vs_pH.txt", "w").write(f"Spearman rho OC vs pH(H2O), LUCAS Romania, n={len(m)}: {rho[0]:.3f} (p={rho[1]:.3f})\n")
print(r.round(3).to_string(index=False)); print(open(out + "/LUCAS_RO_OC_vs_pH.txt").read())
