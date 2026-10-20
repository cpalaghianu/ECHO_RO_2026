"""Step 5 (O4) - how much variation is attributable to the collector versus the place?

A. One-way ICC(1) on mid-ranks (Shrout & Fleiss 1979) with collector, 10-km cell and 25-km cell as
   groups (groups with >= 2 samples); percentile cluster-bootstrap CI (groups resampled).
B. Crossed random-effects model (Gaussian LMM, REML) fitted to the rank-based inverse-normal
   (Blom) transform of each indicator:  y_i = mu + a_collector(i) + b_cell25(i) + c_block100(i) + e_i.
   The 100-km block term absorbs regional gradients, so that they are not attributed to collectors
   (collectors work in limited areas).
   Variance shares of collector, 25-km cell, 100-km block and residual; likelihood-ratio tests of each variance
   component against the model without it, with the 50:50 chi-square(0,1) mixture p-value for a
   variance on the boundary (Self & Liang 1987; Stram & Lee 1994). Collector and place are
   partially confounded (collectors work close to home); the crossed model uses collectors who
   sample in more than one cell to separate them, so shares are descriptive, not causal.
C. Kruskal-Wallis test among collectors with >= 10 samples (effect size epsilon-squared).
"""
import _bootstrap  # noqa: F401
import warnings
import numpy as np, pandas as pd
from scipy import stats
import statsmodels.api as sm
from echo_ro import config as cfg, io  # (not "C": patsy resolves C() in the caller namespace)
from echo_ro.stats_utils import grid_ids, icc1, blom_inverse_normal

warnings.filterwarnings("ignore")
rng = io.rng("m05")
d = io.load()
N = len(d)
cl = pd.factorize(d.collected_by)[0]
cells = {"10-km cell": grid_ids(d.x_km.values, d.y_km.values, 10.0), "25-km cell": grid_ids(d.x_km.values, d.y_km.values, 25.0)}
V = {"pH": d.soil_ph_field.values.astype(float), "Organic matter": d.soil_organic_matter.values.astype(float),
     "Earthworms (counts)": d.earthworms_count.values.astype(float)}

# A. one-way ICC
rows = []
for var, v in V.items():
    ok = ~np.isnan(v); rk = stats.rankdata(v[ok])
    for gname, g in [("collector", cl[ok])] + [(k, c[ok]) for k, c in cells.items()]:
        o = icc1(rk, g); ug = np.unique(g); ix = [np.where(g == u)[0] for u in ug]
        bs = []
        for _ in range(cfg.N_BOOT_KNN):
            s = rng.integers(0, len(ug), len(ug)); ii = np.concatenate([ix[a] for a in s])
            lab = np.concatenate([np.full(len(ix[a]), j) for j, a in enumerate(s)])
            bs.append(icc1(rk[ii], lab))
        gs = pd.Series(g).value_counts()
        rows.append(dict(variable=var, grouping=gname, groups_ge2=int((gs >= 2).sum()), samples_in_groups=int(gs[gs >= 2].sum()),
                         ICC1=o, ci_lo=np.nanquantile(bs, .025), ci_hi=np.nanquantile(bs, .975)))
io.save_table(pd.DataFrame(rows), "T05_icc_oneway.csv")

# B. crossed random effects: collector + 25-km cell + 100-km block (regional scale)
def fit(df, comps):
    vc = {k: f"0 + C({k})" for k in comps}
    m = sm.MixedLM.from_formula("y ~ 1", groups="one", vc_formula=vc, re_formula="0", data=df)
    return m.fit(reml=True, method=["powell"])   # L-BFGS failed to converge for some reduced models
def components(res):
    return dict(zip(list(res.model.exog_vc.names), np.atleast_1d(res.vcomp)))
COMPS = ["block100", "cell25", "collector"]
out = []
for var, v in V.items():
    ok = ~np.isnan(v)
    df = pd.DataFrame({"y": blom_inverse_normal(v[ok]), "collector": cl[ok], "cell25": cells["25-km cell"][ok],
                       "block100": grid_ids(d.x_km.values[ok], d.y_km.values[ok], 100.0), "one": 1})
    full = fit(df, COMPS)
    vcv = components(full)
    tot = sum(vcv.values()) + full.scale
    row = dict(variable=var, n=int(ok.sum()), collectors=int(df.collector.nunique()), cells25=int(df.cell25.nunique()),
               blocks100=int(df.block100.nunique()),
               **{f"share_{k}": vcv[k] / tot for k in COMPS}, share_residual=full.scale / tot)
    for k in COMPS:                                   # LRT of each component, 50:50 chi2(0,1) mixture
        red = fit(df, [c for c in COMPS if c != k])
        lr = max(2 * (full.llf - red.llf), 0)
        row[f"LR_{k}"] = lr
        row[f"p_{k}"] = 0.5 * stats.chi2.sf(lr, 1)
    out.append(row)
    print(var, {k: round(float(x), 3) if isinstance(x, (float, np.floating)) else x for k, x in row.items()}, flush=True)
io.save_table(pd.DataFrame(out), "T05b_variance_partition_crossed.csv")

# C. Kruskal-Wallis among collectors with >= 10 samples
kw = []
vc_ = pd.Series(cl).value_counts()
big = vc_[vc_ >= 10].index
for var, v in V.items():
    groups = [v[(cl == c) & ~np.isnan(v)] for c in big]
    H, p = stats.kruskal(*groups); n = sum(len(g) for g in groups)
    kw.append(dict(variable=var, collectors=len(big), n=n, H=H, p=p, epsilon2=(H - len(big) + 1) / (n - len(big))))
io.save_table(pd.DataFrame(kw), "T05c_kruskal_wallis_collectors.csv")
