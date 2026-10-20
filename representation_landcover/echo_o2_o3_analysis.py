#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
ECHO Romania - O2 (spatial representativeness) and O3/O4b (indicators vs land cover and altitude)

Inputs
  --samples  RO_voronoi_CLC_altitude.xlsx   one row per sample (1,091): indicators, collector, Voronoi area (ha),
                                             CORINE 2018 level-3 code (Code_18), level-1 class (Corine2018), altitude (m)
  --grid     RO_grid_count_area.xlsx         30 x 30 km cells intersecting Romania: No_of_points, area_ha
  --ref-clc  (optional) CSV  class,area_ha   true CORINE level-1 areas of Romania (from the authors' GIS)
  --ref-alt  (optional) CSV  band,area_ha    true area per altitude band (labels = BAND_LABELS below)

Everything that involves inference treats the COLLECTOR as the cluster (cluster bootstrap, mixed models, GEE),
because the collector explains a large share of the variance (O4c).
Pre-specified decisions (fixed before the tests were run) are marked PRE in the comments.
"""
import argparse, json, os, sys, warnings
import numpy as np, pandas as pd
from scipy import stats
import statsmodels.formula.api as smf
import statsmodels.api as sm
from statsmodels.stats.multitest import multipletests

warnings.filterwarnings("ignore")

# ----------------------------------------------------------------------------- constants (PRE)
SEED = 20261019
CLC1 = {1: "Artificial surfaces", 2: "Agricultural areas", 3: "Forest and semi-natural areas",
        4: "Wetlands", 5: "Water bodies"}
TEST_CLC = [1, 2, 3]                     # PRE: tests only for classes with n >= 30; classes 4-5 (n = 11) are descriptive
ALT_EDGES = [-np.inf, 200, 500, 800, 1200, np.inf]          # PRE: altitude bands (m a.s.l.)
BAND_LABELS = ["<200", "200-499", "500-799", "800-1199", ">=1200"]
PH_EDGES = [-np.inf, 6.0, 7.5, np.inf]   # same three pH classes as the LUCAS comparison: <6.0 ; 6.0-7.5 ; >7.5 (right-closed at 7.5)
FOREST3 = {311: "Broad-leaved forest", 312: "Coniferous forest", 313: "Mixed forest"}
IND = ["ph", "om", "worms", "poll"]
IND_LABEL = {"ph": "Soil pH", "om": "Organic matter", "worms": "Earthworm count", "poll": "Visible pollutants (count)"}
MIN_FULL_CELL = 0.50                     # PRE: cell counted as 'substantial' if it has >= 50% of the nominal 90,000 ha
NOMINAL_CELL_HA = 90000.0


# ----------------------------------------------------------------------------- helpers
def load(path_samples, path_grid):
    d = pd.read_excel(path_samples)
    d = d.rename(columns={"sample_id": "id", "soil_ph_fi": "ph", "soil_organ": "om", "earthworms": "worms",
                          "pollutants": "poll", "soil_struc": "structure", "soil_textu": "texture",
                          "collected_": "collector", "voronoi_ar": "vor_ha", "Code_18": "clc3", "Corine2018": "clc1",
                          "altit_1": "alt"})
    g = pd.read_excel(path_grid).rename(columns={"No_of_points": "n", "area_ha": "area_ha"})
    g = g[g["id"].notna()].reset_index(drop=True)      # the file ends with a row of totals (no cell id)
    # --- integrity checks (fail loudly)
    assert d["id"].is_unique, "duplicate sample ids"
    assert len(d) == 1091, f"expected 1091 samples, got {len(d)}"
    assert (d.clc3 // 100 == d.clc1).all(), "CLC level-1 class inconsistent with level-3 code"
    assert d.vor_ha.min() > 0 and d.alt.min() >= 0
    assert d.collector.nunique() == 444
    assert g.n.sum() == len(d), "grid counts do not add up to the number of samples"
    rel = abs(d.vor_ha.sum() - g.area_ha.sum()) / g.area_ha.sum()
    assert rel < 0.001, f"Voronoi areas and grid areas differ by {rel:.4%}"
    d["band"] = pd.cut(d.alt, ALT_EDGES, right=False, labels=BAND_LABELS)
    d["cl"] = d.collector.astype("category").cat.codes.values
    dt = pd.to_datetime(d.sampling_d.astype(str), format="%Y/%m/%d %H:%M:%S.%f")
    per = np.where((dt.dt.year == 2025), "Nov-Dec 2025", np.where(dt.dt.month <= 3, "Jan-Mar 2026",
                   np.where(dt.dt.month == 4, "Apr 2026", "May-Oct 2026")))
    d["period"] = pd.Categorical(per, categories=["Apr 2026", "Nov-Dec 2025", "Jan-Mar 2026", "May-Oct 2026"])
    d["poll_any"] = (d.poll > 0).astype(int)
    d["worm_any"] = (d.worms > 0).astype(int)
    return d, g, rel


def gini(x, w=None):
    x = np.asarray(x, float)
    w = np.ones_like(x) if w is None else np.asarray(w, float)
    o = np.argsort(x)
    x, w = x[o], w[o]
    cw = np.cumsum(w) / w.sum()
    cx = np.cumsum(x * w) / (x * w).sum()
    cw0 = np.r_[0, cw]
    cx0 = np.r_[0, cx]
    return 1 - np.sum((cw0[1:] - cw0[:-1]) * (cx0[1:] + cx0[:-1]))


def wquantile(x, w, q):
    """lower weighted quantile: smallest x whose cumulative weight >= q"""
    o = np.argsort(x, kind="stable")
    x, w = np.asarray(x)[o], np.asarray(w)[o]
    cw = np.cumsum(w) / w.sum()
    return x[np.searchsorted(cw, q - 1e-12)]


def blom(x):
    r = stats.rankdata(x)
    return stats.norm.ppf((r - 3 / 8) / (len(x) + 1 / 4))


def pct(a, lo=2.5, hi=97.5):
    a = np.asarray(a, float)
    a = a[np.isfinite(a)]
    return (np.nan, np.nan) if len(a) < 50 else (np.percentile(a, lo), np.percentile(a, hi))


class ClusterBoot:
    """Resamples collectors with replacement; the same replicates are used for every statistic."""
    def __init__(self, cl, B, seed):
        rng = np.random.default_rng(seed)
        order = np.argsort(cl, kind="stable")
        counts = np.bincount(cl)
        starts = np.r_[0, np.cumsum(counts)]
        G = len(counts)
        self.idx = []
        for _ in range(B):
            g = rng.integers(0, G, G)
            self.idx.append(np.concatenate([order[starts[k]:starts[k + 1]] for k in g]))
        self.B = B

    def ci(self, fn):
        vals = np.array([fn(i) for i in self.idx], float)
        return pct(vals), vals


def eps2_kw(values, groups):
    """Kruskal-Wallis H (tie corrected) and epsilon^2 = (H - k + 1)/(n - k); only groups present are used."""
    gs = [values[groups == g] for g in np.unique(groups)]
    gs = [g for g in gs if len(g) > 0]
    k, n = len(gs), sum(len(g) for g in gs)
    if k < 2 or n <= k:
        return np.nan, np.nan, np.nan, k
    try:
        H, p = stats.kruskal(*gs)
    except ValueError:                      # all values identical
        return np.nan, np.nan, np.nan, k
    return H, (H - k + 1) / (n - k), p, k


def cliff(a, b):
    """Cliff's delta = 2U/(n1 n2) - 1  (P(a>b) - P(a<b))"""
    if len(a) == 0 or len(b) == 0:
        return np.nan
    u = stats.mannwhitneyu(a, b, alternative="two-sided").statistic
    return 2 * u / (len(a) * len(b)) - 1


# ----------------------------------------------------------------------------- A. 30-km grid coverage
def grid_analysis(g, n_samples, rng, nsim=10000):
    out = {}
    g = g.copy()
    g["km2"] = g.area_ha / 100
    g["dens"] = g.n / g.km2 * 100                      # samples per 100 km2 (per cell)
    occ = g.n > 0
    sub = g.area_ha >= MIN_FULL_CELL * NOMINAL_CELL_HA
    out["cells_total"] = int(len(g))
    out["cells_occupied"] = int(occ.sum())
    out["cells_occupied_pct"] = float(100 * occ.mean())
    out["cells_substantial"] = int(sub.sum())
    out["cells_substantial_occupied"] = int((occ & sub).sum())
    out["cells_substantial_occupied_pct"] = float(100 * (occ & sub).sum() / sub.sum())
    out["area_total_km2"] = float(g.km2.sum())
    out["area_in_occupied_cells_pct"] = float(100 * g.km2[occ].sum() / g.km2.sum())
    out["cells_empty"] = int((~occ).sum())
    out["cells_ge2"] = int((g.n >= 2).sum())
    out["cells_ge3"] = int((g.n >= 3).sum())
    out["cells_ge2_pct"] = float(100 * (g.n >= 2).mean())
    o = g[occ]
    out["samples_per_occupied_cell_median"] = float(o.n.median())
    out["samples_per_occupied_cell_iqr"] = [float(o.n.quantile(.25)), float(o.n.quantile(.75))]
    out["samples_per_occupied_cell_max"] = int(o.n.max())
    out["density_per100km2_occupied_median"] = float(o.dens.median())
    out["density_national_per100km2"] = float(100 * g.n.sum() / g.km2.sum())
    # area-weighted Lorenz / Gini of the sample counts over the territory (cells weighted by area)
    out["gini_area_weighted"] = float(gini(g.dens.values, g.km2.values))
    s = g.sort_values("dens", ascending=False)
    cumA = s.km2.cumsum() / s.km2.sum()
    cumN = s.n.cumsum() / s.n.sum()
    out["share_samples_in_densest_10pct_area"] = float(100 * cumN[cumA >= 0.10].iloc[0])
    out["share_samples_in_densest_25pct_area"] = float(100 * cumN[cumA >= 0.25].iloc[0])
    out["share_area_with_zero_samples_pct"] = float(100 * g.km2[~occ].sum() / g.km2.sum())
    # CSR benchmark: the same number of points placed at random over the cells, proportional to area
    p = (g.area_ha / g.area_ha.sum()).values
    sims = rng.multinomial(n_samples, p, size=nsim)
    occ_sim = (sims > 0).sum(1)
    cov_sim = ((sims > 0) * g.km2.values).sum(1) / g.km2.sum() * 100
    exp_occ = float(np.sum(1 - (1 - p) ** n_samples))     # analytic expectation (check against simulation)
    out["csr_expected_occupied_analytic"] = exp_occ
    out["csr_expected_occupied_sim_mean"] = float(occ_sim.mean())
    out["csr_occupied_sim_2_5_97_5"] = [float(np.percentile(occ_sim, 2.5)), float(np.percentile(occ_sim, 97.5))]
    out["csr_min_occupied_sim"] = int(occ_sim.min())
    out["csr_p_obs_le"] = float((np.sum(occ_sim <= occ.sum()) + 1) / (nsim + 1))
    out["csr_ratio_obs_to_expected"] = float(occ.sum() / exp_occ)
    out["csr_area_coverage_sim_mean_pct"] = float(cov_sim.mean())
    out["csr_expected_empty"] = float(len(g) - exp_occ)
    out["csr_nsim"] = nsim
    hist = pd.DataFrame({"samples_in_cell": range(0, 13)})
    hist["observed_cells"] = [(g.n == k).sum() if k < 12 else (g.n >= 12).sum() for k in hist.samples_in_cell]
    hist["csr_expected_cells"] = [np.mean([(sims[i] == k).sum() if k < 12 else (sims[i] >= 12).sum()
                                           for i in range(0, nsim, 20)]) for k in hist.samples_in_cell]
    return out, hist, g


# ----------------------------------------------------------------------------- B. Voronoi footprint
def voronoi_footprint(d, cb):
    a = d.vor_ha.values
    A = a.sum()
    o = {}
    o["total_area_km2"] = float(A / 100)
    o["mean_ha"] = float(a.mean())
    o["median_ha"] = float(np.median(a))
    o["iqr_ha"] = [float(np.percentile(a, 25)), float(np.percentile(a, 75))]
    o["p5_p95_ha"] = [float(np.percentile(a, 5)), float(np.percentile(a, 95))]
    o["max_ha"] = float(a.max())
    o["min_ha"] = float(a.min())
    o["gini_polygon_area"] = float(gini(a))
    o["kish_neff"] = float(A ** 2 / np.sum(a ** 2))
    o["kish_neff_over_n"] = float(o["kish_neff"] / len(a))
    srt = np.sort(a)[::-1]
    k10 = int(np.ceil(0.10 * len(a)))
    o["share_area_in_largest_10pct_polygons"] = float(100 * srt[:k10].sum() / A)
    srt2 = np.sort(a)
    o["share_area_in_smallest_50pct_polygons"] = float(100 * srt2[:len(a) // 2].sum() / A)
    o["n_polygons_lt100ha"] = int((a < 100).sum())
    o["n_polygons_lt1000ha"] = int((a < 1000).sum())
    o["n_polygons_gt50000ha"] = int((a > 50000).sum())
    o["share_area_polygons_gt50000ha_pct"] = float(100 * a[a > 50000].sum() / A)
    ci, _ = cb.ci(lambda i: (a[i].sum() ** 2 / np.sum(a[i] ** 2)) / len(i))
    o["kish_neff_over_n_ci"] = [float(ci[0]), float(ci[1])]
    return o


# ----------------------------------------------------------------------------- C. representation tables
def representation(d, col, levels, labels, cb, ref=None):
    """Shares of samples, Voronoi-estimated shares, relative sampling intensity (RSI) and, if given, reference ratio RR."""
    codes = pd.Categorical(d[col], categories=levels).codes
    ok_all = codes >= 0
    a = d.vor_ha.values
    n_tot = len(d)
    A_tot = a.sum()
    K = len(levels)

    def stats_(idx):
        c = codes[idx]
        m = c >= 0
        nn = np.bincount(c[m], minlength=K).astype(float)
        aa = np.bincount(c[m], weights=a[idx][m], minlength=K)
        return nn, aa, m.sum(), aa.sum()
    nn, aa, nt, At = stats_(np.arange(n_tot))
    rows = []
    boots = [stats_(i) for i in cb.idx]
    for k in range(K):
        share = nn[k] / nt
        vsh = aa[k] / At
        r = {"class": labels[k], "n": int(nn[k]), "share_samples_pct": 100 * share,
             "voronoi_area_km2": aa[k] / 100, "voronoi_share_pct": 100 * vsh,
             "mean_polygon_ha": aa[k] / nn[k] if nn[k] else np.nan,
             "RSI": share / vsh if vsh > 0 else np.nan}
        bs = np.array([b[0][k] / b[2] for b in boots])
        r["share_samples_ci"] = tuple(100 * x for x in pct(bs))
        bv = np.array([(b[1][k] / b[3]) if b[3] > 0 else np.nan for b in boots])
        r["voronoi_share_ci"] = tuple(100 * x for x in pct(bv))
        with np.errstate(divide="ignore", invalid="ignore"):
            br = np.where(bv > 0, bs / bv, np.nan)
        r["RSI_ci"] = pct(br)
        if ref is not None and labels[k] in ref:
            rs = ref[labels[k]] / sum(ref.values())
            r["ref_share_pct"] = 100 * rs
            r["RR_vs_reference"] = share / rs if rs > 0 else np.nan
            r["RR_ci"] = pct(bs / rs)
        rows.append(r)
    return pd.DataFrame(rows)


# ----------------------------------------------------------------------------- D. declustered descriptives
def desc_stats(d, w):
    """Scalar summaries computed with weights w (sum to 1)."""
    o = {}
    ph, om = d.ph.values, d.om.values
    mo = ~np.isnan(om)
    o["pH median"] = wquantile(ph, w, 0.5)
    o["pH Q1"] = wquantile(ph, w, 0.25)
    o["pH Q3"] = wquantile(ph, w, 0.75)
    cls = np.where(ph < 6.0, 0, np.where(ph <= 7.5, 1, 2))
    for k, nm in enumerate(["pH <6.0", "pH 6.0-7.5", "pH >7.5"]):
        o[nm + " (%)"] = 100 * w[cls == k].sum()
    wm = w[mo] / w[mo].sum()
    o["OM median"] = wquantile(om[mo], wm, 0.5)
    for v in [1.5, 2, 2.5, 3.5, 5]:
        o[f"OM class {v} (%)"] = 100 * wm[om[mo] == v].sum()
    o["Earthworms >=1 (%)"] = 100 * w[d.worms.values > 0].sum()
    o["Earthworm count, mean"] = float(np.sum(w * d.worms.values))
    o["Earthworm count, median"] = wquantile(d.worms.values, w, 0.5)
    o["Pollutants >=1 (%)"] = 100 * w[d.poll.values > 0].sum()
    for s in ["friable", "intact", "firm", "compact", "high compact"]:
        o[f"Structure {s} (%)"] = 100 * w[d.structure.values == s].sum()
    mt = d.texture.values != "other"
    wt = w[mt] / w[mt].sum()
    for s in ["sandy loam", "clay loam", "sandy clay", "silty clay", "sandy", "silty", "clayey"]:
        o[f"Texture {s} (%)"] = 100 * wt[d.texture.values[mt] == s].sum()
    return o


def descriptives(d, cb):
    n = len(d)
    w_u = np.full(n, 1 / n)
    w_v = (d.vor_ha / d.vor_ha.sum()).values
    cap = np.percentile(d.vor_ha, 99)
    w_c = np.minimum(d.vor_ha.values, cap)
    w_c = w_c / w_c.sum()                                         # sensitivity: weights capped at the 99th percentile
    base = {k: desc_stats(d, w) for k, w in [("unweighted", w_u), ("voronoi", w_v), ("voronoi_capped99", w_c)]}

    def boot_one(idx, which):
        dd = d.iloc[idx]
        a = dd.vor_ha.values
        if which == "unweighted":
            w = np.full(len(dd), 1 / len(dd))
        elif which == "voronoi":
            w = a / a.sum()
        else:
            w = np.minimum(a, cap) / np.minimum(a, cap).sum()
        return desc_stats(dd, w)
    res = {}
    for which in base:
        reps = [boot_one(i, which) for i in cb.idx[:2000]]        # 2,000 replicates are enough for these CIs
        res[which] = {k: pct([r[k] for r in reps]) for k in base[which]}
    rows = []
    for k in base["unweighted"]:
        row = {"statistic": k}
        for which in base:
            row[which] = base[which][k]
            row[which + "_lo"], row[which + "_hi"] = res[which][k]
        rows.append(row)
    return pd.DataFrame(rows)


def group_descriptives(d):
    rows = []
    for fac, lev, lab in [("clc1", TEST_CLC + [4, 5], CLC1), ("band", BAND_LABELS, None)]:
        for L in lev:
            s = d[d[fac] == L]
            if len(s) == 0:
                continue
            q = lambda x, p: np.nanpercentile(x, p)
            rows.append({"factor": "CORINE level 1" if fac == "clc1" else "Altitude band (m)",
                         "level": CLC1[L] if fac == "clc1" else L, "n": len(s), "collectors": s.collector.nunique(),
                         "pH median": s.ph.median(), "pH IQR": f"{q(s.ph,25):.1f}-{q(s.ph,75):.1f}",
                         "OM median": s.om.median(), "OM IQR": f"{q(s.om,25):.1f}-{q(s.om,75):.1f}",
                         "Earthworms median": s.worms.median(), "Earthworms >=1 (%)": 100 * (s.worms > 0).mean(),
                         "Pollutants >=1 (%)": 100 * (s.poll > 0).mean(),
                         "altitude median (m)": s.alt.median()})
    return pd.DataFrame(rows)


# ----------------------------------------------------------------------------- E. Table 3 statistics
def table3(d, cb):
    vals = {k: d[k].values.astype(float) for k in IND}
    alt = d.alt.values.astype(float)
    # --- E1 Spearman correlations (pairwise complete)
    pairs = [("ph", "om"), ("ph", "worms"), ("ph", "poll"), ("om", "worms"), ("om", "poll"), ("worms", "poll"),
             ("alt", "ph"), ("alt", "om"), ("alt", "worms"), ("alt", "poll")]
    src = dict(vals, alt=alt)

    def rho(idx, a, b):
        x, y = src[a][idx], src[b][idx]
        m = ~(np.isnan(x) | np.isnan(y))
        if m.sum() < 10:
            return np.nan
        return stats.spearmanr(x[m], y[m])[0]
    full = np.arange(len(d))
    rows = []
    for a, b in pairs:
        x, y = src[a], src[b]
        m = ~(np.isnan(x) | np.isnan(y))
        r, p = stats.spearmanr(x[m], y[m])
        ci, _ = cb.ci(lambda i, a=a, b=b: rho(i, a, b))
        rows.append({"test": "Spearman", "x": a, "y": b, "n": int(m.sum()), "stat": r, "ci_lo": ci[0], "ci_hi": ci[1],
                     "p_naive": p, "family": "correlation"})
    # --- E2 Kruskal-Wallis across CLC level 1 (classes 1-3) and altitude bands; epsilon^2 with cluster-bootstrap CI
    for fac, lev in [("clc1", TEST_CLC), ("band", BAND_LABELS)]:
        inlev = d[fac].isin(lev).values
        g = d[fac].astype(str).values
        for k in IND:
            ok = inlev & ~np.isnan(vals[k])
            H, e2, p, kk = eps2_kw(vals[k][ok], g[ok])

            def e2_boot(i, k=k, ok=ok, g=g):
                m = ok[i]
                return eps2_kw(vals[k][i][m], g[i][m])[1]
            ci, _ = cb.ci(e2_boot)
            rows.append({"test": "Kruskal-Wallis", "x": fac, "y": k, "n": int(ok.sum()), "stat": e2, "ci_lo": ci[0],
                         "ci_hi": ci[1], "p_naive": p, "H": H, "k": kk, "family": "kruskal"})
    T = pd.DataFrame(rows)
    T["p_BH_naive"] = multipletests(T.p_naive, method="fdr_bh")[1]        # one family: 10 correlations + 8 KW tests = 18
    return T


def posthoc(d, cb):
    """A-priori contrasts (PRE: forest vs agricultural for pH, OM, earthworms; the user's working expectations)
       and exploratory contrasts. Cliff's delta (first minus second group) with cluster-bootstrap CI."""
    contrasts = [
        ("Forest/semi-natural vs agricultural", "clc1", 3, "clc1", 2, "PRE"),
        ("Forest/semi-natural vs artificial", "clc1", 3, "clc1", 1, "exploratory"),
        ("Agricultural vs artificial", "clc1", 2, "clc1", 1, "exploratory"),
        ("Coniferous vs broad-leaved forest", "clc3", 312, "clc3", 311, "exploratory"),
        ("Altitude >=800 m vs <200 m", "alt_hi", 1, "alt_lo", 1, "exploratory"),
    ]
    d = d.copy()
    d["alt_hi"] = (d.alt >= 800).astype(int)
    d["alt_lo"] = (d.alt < 200).astype(int)
    rows = []
    for name, f1, v1, f2, v2, kind in contrasts:
        m1 = (d[f1] == v1).values
        m2 = (d[f2] == v2).values
        for k in IND:
            y = d[k].values.astype(float)

            def dl(i, y=y, m1=m1, m2=m2):
                a = y[i][m1[i]]
                b = y[i][m2[i]]
                return cliff(a[~np.isnan(a)], b[~np.isnan(b)])
            a = y[m1]; a = a[~np.isnan(a)]
            b = y[m2]; b = b[~np.isnan(b)]
            p = stats.mannwhitneyu(a, b, alternative="two-sided").pvalue
            ci, _ = cb.ci(dl)
            rows.append({"contrast": name, "status": kind, "indicator": k, "n1": len(a), "n2": len(b),
                         "median1": np.median(a), "median2": np.median(b), "cliffs_delta": cliff(a, b),
                         "ci_lo": ci[0], "ci_hi": ci[1], "p_naive": p})
    R = pd.DataFrame(rows)
    R["p_BH_naive"] = multipletests(R.p_naive, method="fdr_bh")[1]
    return R


def adjusted_models(d):
    """Collector random intercept + CORINE (3 classes, reference agricultural) + altitude (per 100 m).
       pH, OM, earthworms: Blom normal scores, LMM (ML for LRT, REML for estimates). Pollutants: presence, logistic GEE."""
    dd = d[d.clc1.isin(TEST_CLC)].copy()
    dd["alt100"] = dd.alt / 100
    dd["clc"] = pd.Categorical(dd.clc1, categories=[2, 1, 3])
    rows, icc = [], {}
    for k in ["ph", "om", "worms"]:
        x = dd.dropna(subset=[k]).copy()
        x["y"] = blom(x[k].values)
        full_ml = smf.mixedlm("y ~ C(clc) + alt100", x, groups=x.collector).fit(reml=False)
        no_clc = smf.mixedlm("y ~ alt100", x, groups=x.collector).fit(reml=False)
        no_alt = smf.mixedlm("y ~ C(clc)", x, groups=x.collector).fit(reml=False)
        null = smf.mixedlm("y ~ 1", x, groups=x.collector).fit(reml=True)
        icc[k] = float(null.cov_re.iloc[0, 0] / (null.cov_re.iloc[0, 0] + null.scale))
        m = smf.mixedlm("y ~ C(clc) + alt100", x, groups=x.collector).fit(reml=True)
        lr_clc = 2 * (full_ml.llf - no_clc.llf)
        lr_alt = 2 * (full_ml.llf - no_alt.llf)
        ci = m.conf_int()
        for term in m.fe_params.index:
            if term == "Intercept":
                continue
            lab = {"C(clc)[T.1]": "Artificial vs agricultural", "C(clc)[T.3]": "Forest/semi-natural vs agricultural",
                   "alt100": "Altitude, per +100 m"}[term]
            rows.append({"indicator": k, "term": lab, "n": len(x), "collectors": x.collector.nunique(),
                         "estimate": m.fe_params[term], "ci_lo": ci.loc[term, 0], "ci_hi": ci.loc[term, 1],
                         "p_wald": m.pvalues[term], "model": "LMM, Blom score (SD units)",
                         "LRT_clc_chi2_df2": lr_clc, "LRT_clc_p": stats.chi2.sf(lr_clc, 2),
                         "LRT_alt_chi2_df1": lr_alt, "LRT_alt_p": stats.chi2.sf(lr_alt, 1)})
    x = dd.copy()
    gee = smf.gee("poll_any ~ C(clc) + alt100", "collector", x, family=sm.families.Binomial(),
                  cov_struct=sm.cov_struct.Exchangeable()).fit()
    ci = gee.conf_int()
    for term in gee.params.index:
        if term == "Intercept":
            continue
        lab = {"C(clc)[T.1]": "Artificial vs agricultural", "C(clc)[T.3]": "Forest/semi-natural vs agricultural",
               "alt100": "Altitude, per +100 m"}[term]
        rows.append({"indicator": "poll", "term": lab, "n": len(x), "collectors": x.collector.nunique(),
                     "estimate": np.exp(gee.params[term]), "ci_lo": np.exp(ci.loc[term, 0]),
                     "ci_hi": np.exp(ci.loc[term, 1]), "p_wald": gee.pvalues[term],
                     "model": "Logistic GEE, collector clusters (odds ratio)"})
    R = pd.DataFrame(rows)
    R["p_BH"] = multipletests(R.p_wald, method="fdr_bh")[1]
    # unadjusted collector-level ICC (single random intercept; the crossed model of O4c is the reference estimate)
    return R, icc



def sensitivity_models(d):
    """Sensitivity of the adjusted models: (S1) + sampling period (confounding of altitude with the season of sampling),
       (S2) + altitude^2 (non-linearity; LRT), for the three Blom-score indicators."""
    dd = d[d.clc1.isin(TEST_CLC)].copy()
    dd["alt100"] = dd.alt / 100
    dd["alt100_2"] = dd.alt100 ** 2
    dd["clc"] = pd.Categorical(dd.clc1, categories=[2, 1, 3])
    rows = []
    rho_alt_day = stats.spearmanr(dd.alt, (pd.to_datetime(dd.sampling_d.astype(str), format="%Y/%m/%d %H:%M:%S.%f")
                                           - pd.Timestamp("2025-11-01")).dt.days)
    for k in ["ph", "om", "worms"]:
        x = dd.dropna(subset=[k]).copy()
        x["y"] = blom(x[k].values)
        base = smf.mixedlm("y ~ C(clc) + alt100", x, groups=x.collector).fit(reml=False)
        per = smf.mixedlm("y ~ C(clc) + alt100 + C(period)", x, groups=x.collector).fit(reml=False)
        quad = smf.mixedlm("y ~ C(clc) + alt100 + alt100_2", x, groups=x.collector).fit(reml=False)
        per_r = smf.mixedlm("y ~ C(clc) + alt100 + C(period)", x, groups=x.collector).fit(reml=True)
        ci = per_r.conf_int()
        for term, lab in [("C(clc)[T.1]", "Artificial vs agricultural"), ("C(clc)[T.3]", "Forest/semi-natural vs agricultural"),
                          ("alt100", "Altitude, per +100 m")]:
            rows.append({"indicator": k, "model": "S1: + sampling period (REML)", "term": lab, "estimate": per_r.fe_params[term],
                         "ci_lo": ci.loc[term, 0], "ci_hi": ci.loc[term, 1], "p_wald": per_r.pvalues[term],
                         "LRT_period_chi2_df3": 2 * (per.llf - base.llf), "LRT_period_p": stats.chi2.sf(2 * (per.llf - base.llf), 3),
                         "LRT_alt2_chi2_df1": np.nan, "LRT_alt2_p": np.nan})
        rows.append({"indicator": k, "model": "S2: + altitude^2 (ML)", "term": "Altitude^2 (per 100 m, squared)",
                     "estimate": quad.fe_params["alt100_2"], "ci_lo": np.nan, "ci_hi": np.nan, "p_wald": quad.pvalues["alt100_2"],
                     "LRT_period_chi2_df3": np.nan, "LRT_period_p": np.nan,
                     "LRT_alt2_chi2_df1": 2 * (quad.llf - base.llf), "LRT_alt2_p": stats.chi2.sf(2 * (quad.llf - base.llf), 1)})
    return pd.DataFrame(rows), {"spearman_altitude_vs_day_of_campaign": [float(rho_alt_day[0]), float(rho_alt_day[1])],
                                "period_counts": d.period.value_counts().to_dict(),
                                "altitude_median_by_period": d.groupby("period").alt.median().to_dict()}

# ----------------------------------------------------------------------------- figures
def figures(outdir, reprC, reprA, hist, d):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    plt.rcParams.update({"font.family": "DejaVu Sans", "font.size": 8, "axes.spines.top": False,
                         "axes.spines.right": False, "axes.edgecolor": "#52514e", "axes.labelcolor": "#0b0b0b",
                         "xtick.color": "#52514e", "ytick.color": "#52514e"})
    B1, B2, B3, INK, MUT = "#2a78d6", "#eb6834", "#1baf7a", "#0b0b0b", "#52514e"
    fig, ax = plt.subplots(1, 3, figsize=(10.8, 3.5), gridspec_kw={"width_ratios": [1.05, 1.25, 1.25]})

    def bars(a, df, title, short):
        x = np.arange(len(df))
        has_ref = "ref_share_pct" in df
        w = 0.27 if has_ref else 0.36
        o1, o2 = (-w, 0.0) if has_ref else (-w / 2, w / 2)
        a.bar(x + o1, df.share_samples_pct, w, color=B1, label="Share of samples", zorder=3)
        a.bar(x + o2, df.voronoi_share_pct, w, color=B2, label="Share of Voronoi area", zorder=3)
        if has_ref:
            a.bar(x + w, df.ref_share_pct, w, color=B3, label="Share of true area", zorder=3)
        for i, r in enumerate(df.itertuples()):
            a.errorbar(i + o1, r.share_samples_pct, yerr=[[r.share_samples_pct - r.share_samples_ci[0]],
                                                          [r.share_samples_ci[1] - r.share_samples_pct]],
                       color=INK, lw=0.8, capsize=2, zorder=4)
            a.text(i + o1, r.share_samples_ci[1] + 1.2, f"{r.share_samples_pct:.0f}", ha="center", fontsize=7, color=INK)
            a.text(i + o2, r.voronoi_share_pct + 1.2, f"{r.voronoi_share_pct:.0f}", ha="center", fontsize=7, color=INK)
            if has_ref:
                a.text(i + w, r.ref_share_pct + 1.2, f"{r.ref_share_pct:.0f}", ha="center", fontsize=7, color=INK)
        a.set_xticks(x)
        a.set_xticklabels(short)
        a.set_ylabel("%")
        a.set_title(title, loc="left", fontsize=8.5, color=INK)
        a.grid(axis="y", color="#e6e5e1", lw=0.6, zorder=0)
    c = reprC[reprC["class"].isin([CLC1[1], CLC1[2], CLC1[3], CLC1[4], CLC1[5]])]
    bars(ax[0], c, "(a) CORINE level 1", ["Artif.", "Agric.", "Forest/\nsemi-nat.", "Wetl.\n(n=2)", "Water\n(n=9)"])
    ax[0].set_ylim(0, 72)
    ax[0].legend(frameon=False, fontsize=6.5, loc="upper right")
    bars(ax[1], reprA, "(b) Altitude band (m a.s.l.)", [f"{b}\n(n={n})" for b, n in zip(reprA["class"], reprA.n)])
    ax[1].legend(frameon=False, fontsize=7, loc="upper right")
    h = hist
    x = np.arange(len(h))
    ax[2].bar(x - 0.2, h.observed_cells, 0.4, color=B1, label="Observed", zorder=3)
    ax[2].bar(x + 0.2, h.csr_expected_cells, 0.4, color=MUT, alpha=0.55, label="Random placement", zorder=3)
    ax[2].set_xticks(x)
    ax[2].set_xticklabels([str(i) if i < 12 else "12+" for i in h.samples_in_cell], fontsize=7)
    ax[2].set_xlabel("Samples per 30 × 30 km cell")
    ax[2].set_ylabel("Number of cells")
    ax[2].set_title("(c) Cell occupancy", loc="left", fontsize=8.5, color=INK)
    ax[2].legend(frameon=False, fontsize=7)
    ax[2].grid(axis="y", color="#e6e5e1", lw=0.6, zorder=0)
    fig.tight_layout()
    fig.savefig(os.path.join(outdir, "Fig_O2_representation.png"), dpi=300)
    fig.savefig(os.path.join(outdir, "Fig_O2_representation.pdf"))
    plt.close(fig)

    # O3: indicators by CORINE class and altitude band
    fig, ax = plt.subplots(2, 3, figsize=(10.2, 5.6))
    sets = [("clc1", [1, 2, 3], ["Artificial", "Agricultural", "Forest/\nsemi-nat."]),
            ("band", BAND_LABELS, [f"{b}" for b in BAND_LABELS])]
    cols = {"ph": "Soil pH", "om": "Organic matter class", "worms": "Earthworm count"}
    for r, (fac, lev, lab) in enumerate(sets):
        for c_, k in enumerate(["ph", "om", "worms"]):
            a = ax[r, c_]
            data = [d.loc[d[fac] == L, k].dropna().values for L in lev]
            bp = a.boxplot(data, widths=0.55, patch_artist=True, showfliers=(k == "worms"),
                           medianprops=dict(color=INK, lw=1.2), whiskerprops=dict(color=MUT, lw=0.8),
                           capprops=dict(color=MUT, lw=0.8), flierprops=dict(marker="o", ms=2, mfc="none", mec=MUT))
            for patch in bp["boxes"]:
                patch.set_facecolor("#cfe0f6")
                patch.set_edgecolor(B1)
            a.set_xticks(range(1, len(lev) + 1))
            a.set_xticklabels([f"{l}\n(n={len(x)})" for l, x in zip(lab, data)], fontsize=7)
            a.set_ylabel(cols[k])
            if k == "worms":
                a.set_ylim(-0.5, 16)
            a.grid(axis="y", color="#e6e5e1", lw=0.6, zorder=0)
            if r == 0:
                a.set_title(f"({'abc'[c_]}) {cols[k]} by CORINE class", loc="left", fontsize=8.5, color=INK)
            else:
                a.set_title(f"({'def'[c_]}) {cols[k]} by altitude (m)", loc="left", fontsize=8.5, color=INK)
    fig.tight_layout()
    fig.savefig(os.path.join(outdir, "Fig_O3_indicators_by_class_altitude.png"), dpi=300)
    fig.savefig(os.path.join(outdir, "Fig_O3_indicators_by_class_altitude.pdf"))
    plt.close(fig)


# ----------------------------------------------------------------------------- main
def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--samples", required=True)
    ap.add_argument("--grid", required=True)
    ap.add_argument("--out", required=True)
    ap.add_argument("--B", type=int, default=4000)
    ap.add_argument("--ref-clc")
    ap.add_argument("--ref-alt")
    a = ap.parse_args()
    os.makedirs(a.out, exist_ok=True)
    rng = np.random.default_rng(SEED)
    d, g, rel = load(a.samples, a.grid)
    cb = ClusterBoot(d.cl.values, a.B, SEED + 1)
    S = {"seed": SEED, "B": a.B, "n_samples": len(d), "n_collectors": int(d.collector.nunique()),
         "voronoi_vs_grid_area_rel_diff": float(rel)}

    grid_s, hist, gdf = grid_analysis(g, len(d), rng)
    S["grid"] = grid_s
    hist.to_csv(os.path.join(a.out, "O2_grid_occupancy_hist.csv"), index=False)
    S["voronoi"] = voronoi_footprint(d, cb)

    refC = refA = None
    if a.ref_clc:
        r = pd.read_csv(a.ref_clc)
        refC = {CLC1[int(k)] if str(k).isdigit() else k: v for k, v in zip(r.iloc[:, 0], r.iloc[:, 1])}
    if a.ref_alt:
        r = pd.read_csv(a.ref_alt)
        refA = {str(k): v for k, v in zip(r.iloc[:, 0], r.iloc[:, 1])}
    reprC = representation(d, "clc1", [1, 2, 3, 4, 5], [CLC1[i] for i in range(1, 6)], cb, refC)
    reprA = representation(d, "band", BAND_LABELS, BAND_LABELS, cb, refA)
    l3 = sorted(d.clc3.unique())
    reprL3 = representation(d, "clc3", l3, [str(x) for x in l3], cb)
    def fmt(df):
        df = df.copy()
        for c in df.columns:
            if df[c].map(lambda x: isinstance(x, tuple)).any():
                df[c] = df[c].map(lambda x: f"{x[0]:.2f} to {x[1]:.2f}" if isinstance(x, tuple) else x)
        return df
    for nm, df in [("O2_representation_clc1", reprC), ("O2_representation_altitude", reprA),
                   ("O2_representation_clc3_supplementary", reprL3)]:
        fmt(df).to_csv(os.path.join(a.out, nm + ".csv"), index=False)

    desc = descriptives(d, cb)
    desc.to_csv(os.path.join(a.out, "O3_descriptive_national_unweighted_vs_voronoi.csv"), index=False)
    gd = group_descriptives(d)
    gd.to_csv(os.path.join(a.out, "O3_descriptive_by_group.csv"), index=False)
    T = table3(d, cb)
    T.to_csv(os.path.join(a.out, "O3_table3_tests.csv"), index=False)
    P = posthoc(d, cb)
    P.to_csv(os.path.join(a.out, "O3_posthoc_contrasts.csv"), index=False)
    M, icc = adjusted_models(d)
    M.to_csv(os.path.join(a.out, "O3_adjusted_models.csv"), index=False)
    S["icc_collector_null_models"] = icc
    SM, sinfo = sensitivity_models(d)
    SM.to_csv(os.path.join(a.out, "O3_sensitivity_models.csv"), index=False)
    S["season_vs_altitude"] = sinfo
    # altitude by class (explains the confounding)
    S["altitude_median_by_clc"] = d.groupby("clc1").alt.median().to_dict()
    S["share_forest_by_band_pct"] = (100 * d.groupby("band")["clc1"].apply(lambda s: (s == 3).mean())).to_dict()
    S["max_altitude_m"] = int(d.alt.max())
    S["n_ge_1500m"] = int((d.alt >= 1500).sum())
    S["collectors_multiclass"] = int((d.groupby("collector").clc1.nunique() > 1).sum())
    S["collectors_multi_sample"] = int((d.groupby("collector").size() > 1).sum())
    json.dump(S, open(os.path.join(a.out, "summary.json"), "w"), indent=1, default=float)
    figures(a.out, reprC, reprA, hist, d)

    with pd.ExcelWriter(os.path.join(a.out, "ECHO_RO_O2_O3_results.xlsx")) as xw:
        pd.DataFrame([{"key": k, "value": json.dumps(v) if isinstance(v, (dict, list)) else v}
                      for k, v in S.items() if k not in ("grid", "voronoi")] +
                     [{"key": "grid." + k, "value": json.dumps(v) if isinstance(v, list) else v} for k, v in S["grid"].items()] +
                     [{"key": "voronoi." + k, "value": json.dumps(v) if isinstance(v, list) else v} for k, v in S["voronoi"].items()]
                     ).to_excel(xw, sheet_name="summary", index=False)
        for nm, df in [("repr_clc1", reprC), ("repr_altitude", reprA), ("repr_clc3", reprL3), ("desc_national", desc),
                       ("desc_groups", gd), ("table3_tests", T), ("posthoc", P), ("adjusted_models", M),
                       ("sensitivity_models", SM), ("grid_hist", hist)]:
            df2 = df.copy()
            for c in df2.columns:
                if df2[c].map(lambda x: isinstance(x, tuple)).any():
                    df2[c] = df2[c].map(lambda x: f"{x[0]:.2f} to {x[1]:.2f}" if isinstance(x, tuple) else x)
            df2.to_excel(xw, sheet_name=nm, index=False)
    print("done ->", a.out)


if __name__ == "__main__":
    main()
