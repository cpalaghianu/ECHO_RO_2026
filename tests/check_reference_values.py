"""Compare key outputs with the reference values of the published run (tolerances allow for
Monte-Carlo noise when ECHO_FAST=1 or when software versions differ). Run after run_all.py."""
import sys, json
from pathlib import Path
import pandas as pd
ROOT = Path(__file__).resolve().parents[1]; T = ROOT / "outputs" / "tables"
ok = True
def check(name, value, ref, tol):
    global ok
    good = abs(value - ref) <= tol
    ok &= good
    print(("PASS" if good else "FAIL"), f"{name}: {value:.4f} (reference {ref}, tol {tol})")
R1 = json.load(open(T / "T01_contribution_summary.json")); R2 = json.load(open(T / "T02_quality_summary.json"))
check("Gini", R1["gini"]["value"], 0.4540, 1e-3)
check("single-sample collectors %", R1["single"]["pct"], 52.25, 0.05)
check("pH heaping ratio", R2["ph_heaping_ratio"]["value"], 1.658, 2e-3)
check("earthworm heaping at 5 (ratio)", R2["earthworms"]["heaping_5"]["ratio"], 2.150, 5e-3)
sp = pd.read_csv(T / "T03_spatial_autocorrelation.csv")
g = sp[(sp.variable == "pH") & (sp["null"] == "global") & (sp.pairs == "all") & (sp.band == "2-5 km")].statistic.iloc[0]
check("Moran I pH 2-5 km (all pairs)", g, 0.3825, 1e-3)
kn = pd.read_csv(T / "T04_knn_consistency.csv")
r = kn[(kn.variable == "pH") & (kn["mode"] == "corrected") & (kn.d0_km == 2.0) & (kn.k == 5)].rho.iloc[0]
check("kNN rho pH (corrected, k=5, d0=2)", r, 0.271, 2e-3)
vp = pd.read_csv(T / "T05b_variance_partition_crossed.csv").set_index("variable")
check("collector variance share, pH", vp.loc["pH", "share_collector"], 0.474, 5e-3)
sys.exit(0 if ok else 1)
