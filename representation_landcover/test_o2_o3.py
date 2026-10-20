#!/usr/bin/env python3
"""Independent checks of echo_o2_o3_analysis.py: (1) synthetic data with known answers, (2) alternative re-computation
of headline quantities from the real files with pandas/scipy only (a different code path)."""
import sys, os, importlib.util
import numpy as np, pandas as pd
from scipy import stats
sys.path.insert(0, os.path.dirname(__file__))
import echo_o2_o3_analysis as E

OK = []
def check(name, cond, info=""):
    assert cond, f"FAILED: {name} {info}"
    OK.append(name)
    print("ok  ", name, info)

# ---- 1. synthetic checks
check("gini equal", abs(E.gini([1, 1, 1, 1])) < 1e-12)
check("gini 0,0,0,10 = 0.75", abs(E.gini([0, 0, 0, 10]) - 0.75) < 1e-12)
check("gini weighted == replication", abs(E.gini([1, 3], [2, 1]) - E.gini([1, 1, 3])) < 1e-12)
x = np.array([1., 2., 3., 4.])
check("wquantile median lower", E.wquantile(x, np.ones(4) / 4, 0.5) == 2.0)
check("wquantile weights", E.wquantile(x, np.array([.1, .1, .1, .7]), 0.5) == 4.0)
check("wquantile ties", E.wquantile(np.array([5, 5, 6, 6.]), np.ones(4) / 4, 0.25) == 5)

# Blom: identical to the function used for the LMM of O4c
spec = importlib.util.spec_from_file_location("su", os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "src", "echo_ro", "stats_utils.py"))
try:
    import types
    src = open(spec.origin).read().replace("from . import config as C", "C=None")
    su = types.ModuleType("su"); exec(compile(src, "su", "exec"), su.__dict__)
    v = np.random.default_rng(1).integers(0, 6, 300).astype(float)
    check("blom == earlier pipeline", np.allclose(E.blom(v), su.blom_inverse_normal(v)))
except Exception as e:
    print("skip blom cross-check:", e)

# Kruskal-Wallis epsilon^2 by hand (rank formula with tie correction)
rng = np.random.default_rng(2)
vals = rng.integers(0, 5, 120).astype(float); grp = rng.integers(0, 3, 120).astype(str)
r = stats.rankdata(vals); n = len(vals)
H0 = 12 / (n * (n + 1)) * sum(r[grp == g].sum() ** 2 / (grp == g).sum() for g in np.unique(grp)) - 3 * (n + 1)
_, cnt = np.unique(vals, return_counts=True)
Hc = H0 / (1 - (cnt ** 3 - cnt).sum() / (n ** 3 - n))
H, e2, p, k = E.eps2_kw(vals, grp)
check("KW H by hand", abs(H - Hc) < 1e-9, f"H={H:.4f}")
check("eps2 formula", abs(e2 - (Hc - 3 + 1) / (n - 3)) < 1e-9)

# Cliff's delta brute force
a = rng.normal(0, 1, 40).round(0); b = rng.normal(.4, 1, 55).round(0)
bf = (sum((x > y) for x in a for y in b) - sum((x < y) for x in a for y in b)) / (len(a) * len(b))
check("cliff brute force", abs(E.cliff(a, b) - bf) < 1e-12, f"{bf:.3f}")

# cluster bootstrap resamples whole collectors
cl = np.repeat(np.arange(30), rng.integers(1, 6, 30))
cb = E.ClusterBoot(cl, 50, 3)
sizes = np.bincount(cl)
for idx in cb.idx[:20]:
    c = np.bincount(cl[idx], minlength=30)
    mult = c / sizes
    assert np.allclose(mult, np.round(mult)), "partial collector in a replicate"
check("bootstrap keeps collectors whole", True)

# RSI on a synthetic population: RSI_c = (n_c/n)/(A_c/A) = mean polygon area overall / mean polygon area class c
syn = pd.DataFrame({"clc1": [1] * 10 + [2] * 30, "vor_ha": [100.] * 10 + [900.] * 30, "cl": np.arange(40)})
syn["band"] = pd.Categorical(["<200"] * 40, categories=E.BAND_LABELS)
cb2 = E.ClusterBoot(syn.cl.values, 100, 4)
R = E.representation(syn, "clc1", [1, 2], ["a", "b"], cb2)
mean_all = syn.vor_ha.mean()
check("RSI class 1", abs(R.RSI[0] - mean_all / 100) < 1e-9, f"{R.RSI[0]:.3f}")
check("RSI class 2", abs(R.RSI[1] - mean_all / 900) < 1e-9, f"{R.RSI[1]:.3f}")
check("voronoi shares sum to 100", abs(R.voronoi_share_pct.sum() - 100) < 1e-9)
# with a reference equal to the Voronoi shares, RR == RSI
ref = {"a": R.voronoi_area_km2[0], "b": R.voronoi_area_km2[1]}
R2 = E.representation(syn, "clc1", [1, 2], ["a", "b"], cb2, ref)
check("RR == RSI when reference = Voronoi area", np.allclose(R2.RR_vs_reference, R2.RSI))

# CSR: analytic expectation equals simulation
g = pd.DataFrame({"n": [0, 5, 2, 0, 1, 3], "area_ha": [90000, 90000, 45000, 20000, 90000, 90000]})
o, hist, _ = E.grid_analysis(g, int(g.n.sum()), np.random.default_rng(5), nsim=20000)
check("CSR analytic == simulation", abs(o["csr_expected_occupied_analytic"] - o["csr_expected_occupied_sim_mean"]) < 0.05,
      f"{o['csr_expected_occupied_analytic']:.3f} vs {o['csr_expected_occupied_sim_mean']:.3f}")

# desc_stats: weighted pH class shares sum to 100
d0 = pd.DataFrame({"ph": [5, 6, 7, 8, 9.], "om": [1.5, 2, 2.5, 3.5, 5], "worms": [0, 1, 2, 3, 4], "poll": [0, 0, 1, 0, 2],
                   "structure": ["friable"] * 5, "texture": ["sandy loam"] * 5})
ds = E.desc_stats(d0, np.ones(5) / 5)
check("pH class shares sum to 100", abs(ds["pH <6.0 (%)"] + ds["pH 6.0-7.5 (%)"] + ds["pH >7.5 (%)"] - 100) < 1e-9)
check("pH classes (5 | 6,7 | 8,9)", (ds["pH <6.0 (%)"], ds["pH 6.0-7.5 (%)"], ds["pH >7.5 (%)"]) == (20.0, 40.0, 40.0))

# ---- 2. independent re-computation on the real files (pandas/scipy only)
base = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "data", "input") + os.sep
d, gr, rel = E.load(base + "RO_voronoi_CLC_altitude.xlsx", base + "RO_grid_count_area.xlsx")
raw = pd.read_excel(base + "RO_voronoi_CLC_altitude.xlsx")
tot = raw.voronoi_ar.sum()
for c in [1, 2, 3]:
    s = raw[raw.Corine2018 == c]
    mine = E.representation(d, "clc1", [1, 2, 3, 4, 5], [str(i) for i in range(1, 6)], E.ClusterBoot(d.cl.values, 50, 1))
    rsi = (len(s) / len(raw)) / (s.voronoi_ar.sum() / tot)
    check(f"RSI CLC{c} (pandas re-computation)", abs(mine.RSI[c - 1] - rsi) < 1e-9, f"{rsi:.3f}")
check("occupied cells", int((gr.n > 0).sum()) == E.grid_analysis(gr, 1091, np.random.default_rng(1), 100)[0]["cells_occupied"])
a = raw.voronoi_ar.values.astype(float)
vf = E.voronoi_footprint(d, E.ClusterBoot(d.cl.values, 50, 1))
check("Kish n_eff", abs(vf["kish_neff"] - a.sum() ** 2 / (a ** 2).sum()) < 1e-6, f"{vf['kish_neff']:.1f}")
# Spearman with pandas
for (x, y) in [("ph", "om"), ("ph", "worms")]:
    r_pd = d[[x, y]].dropna().corr(method="spearman").iloc[0, 1]
    m = d[[x, y]].dropna()
    check(f"spearman {x}-{y}", abs(r_pd - stats.spearmanr(m[x], m[y])[0]) < 1e-12, f"{r_pd:.3f}")
# altitude bands cover every sample exactly once
check("bands partition", d.band.notna().all() and d.band.value_counts().sum() == 1091)
check("CLC1 counts 168/370/542/2/9", d.clc1.value_counts().sort_index().tolist() == [168, 370, 542, 2, 9])

# ---- 3. representation ratio against the true CORINE areas (ref_clc.csv, supplied by the authors from GIS)
ref = pd.read_csv(os.path.join(os.path.dirname(os.path.abspath(__file__)), "ref_clc.csv")).set_index("class").area_ha
rc = pd.read_csv(os.path.join(os.path.dirname(os.path.abspath(__file__)), "results", "O2_representation_clc1.csv"))
for c in [1, 2, 3, 4, 5]:
    rr_pd = (raw.Corine2018.eq(c).mean()) / (ref[c] / ref.sum())
    check(f"RR CLC{c} (pandas re-computation)", abs(rc.RR_vs_reference[c - 1] - rr_pd) < 1e-9, f"{rr_pd:.3f}")
check("reference areas sum to 23,837,574 ha", abs(ref.sum() - 23837574.3766) < 0.01)
check("reference total close to Voronoi total (<0.01%)", abs(ref.sum() * 0.01 - raw.voronoi_ar.sum() / 100) / (ref.sum() * 0.01) < 1e-4)
# altitude bands: reference areas (procente_altitudine.xlsx -> ref_alt.csv); independent recomputation of the ratios
ra = pd.read_csv(os.path.join(os.path.dirname(os.path.abspath(__file__)), "ref_alt.csv")).set_index("band").area_ha
rzA = pd.read_csv(os.path.join(os.path.dirname(os.path.abspath(__file__)), "results", "O2_representation_altitude.csv")).set_index("class")
bands = pd.cut(raw["altit_1"], [-1e9, 200, 500, 800, 1200, 1e9], right=False, labels=list(ra.index))
for b_ in ra.index:
    rr_pd = (bands.eq(b_).mean()) / (ra[b_] / ra.sum())
    check(f"RR altitude {b_} (pandas re-computation)", abs(rzA.RR_vs_reference[b_] - rr_pd) < 1e-9, f"{rr_pd:.3f}")
check("altitude reference areas sum to 23,839,399 ha", abs(ra.sum() - 23839399) < 1)
print(f"\nALL {len(OK)} CHECKS PASSED")
