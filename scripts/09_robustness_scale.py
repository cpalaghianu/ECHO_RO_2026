"""Step 9 (B2) - are the collector-attributable variance shares an artefact of the rank-based (Blom) normalisation?

The crossed Gaussian LMM of step 5 (collector + 25-km cell + 100-km block + residual, REML) is refitted to the
same indicators on other scales: the original values, mid-ranks rescaled to z-scores (ranks without
normalisation) and, for earthworm counts, log(1 + count) and sqrt(count). The shares are scale-free ratios of
variance components, so they can be compared directly across scales. Also reported: one-way ICC(1) of the
collector on the original values, and a latent-scale check (the share is recomputed on the probit-type
'ordered-category midpoint' score of the ordinal indicators = the Blom score of the category *midranks*, which is
what the reference analysis already uses; hence only the raw scale and the rank scale are alternatives here).
Output: outputs/tables/T09_variance_share_by_scale.csv
"""
import _bootstrap  # noqa: F401
import warnings
import numpy as np, pandas as pd
from scipy import stats
import statsmodels.api as sm
from echo_ro import config as cfg, io
from echo_ro.stats_utils import grid_ids, blom_inverse_normal, icc1

warnings.filterwarnings("ignore")
d = io.load()
cl = pd.factorize(d.collected_by)[0]
cell25 = grid_ids(d.x_km.values, d.y_km.values, 25.0)
block100 = grid_ids(d.x_km.values, d.y_km.values, 100.0)
COMPS = ["block100", "cell25", "collector"]

def fit(df, comps):
    vc = {k: f"0 + C({k})" for k in comps}
    m = sm.MixedLM.from_formula("y ~ 1", groups="one", vc_formula=vc, re_formula="0", data=df)
    return m.fit(reml=True, method=["powell"])

def shares(y, ok):
    df = pd.DataFrame({"y": y[ok], "collector": cl[ok], "cell25": cell25[ok], "block100": block100[ok], "one": 1})
    full = fit(df, COMPS)
    vcv = dict(zip(list(full.model.exog_vc.names), np.atleast_1d(full.vcomp)))
    tot = sum(vcv.values()) + full.scale
    red = fit(df, ["block100", "cell25"])
    lr = max(2 * (full.llf - red.llf), 0)
    return dict(n=int(ok.sum()), share_collector=vcv["collector"] / tot, share_cell25=vcv["cell25"] / tot,
                share_block100=vcv["block100"] / tot, share_residual=full.scale / tot,
                LR_collector=lr, p_collector=0.5 * stats.chi2.sf(lr, 1))

def midrank_z(v):
    r = stats.rankdata(v); z = (r - r.mean()) / r.std(); return z

IND = {"pH": d.soil_ph_field.values.astype(float),
       "Organic matter": d.soil_organic_matter.values.astype(float),
       "Earthworms (counts)": d.earthworms_count.values.astype(float)}
rows = []
for var, v in IND.items():
    ok = ~np.isnan(v)
    scales = {"Blom normal scores (reference)": lambda x: blom_inverse_normal(x),
              "mid-ranks as z-scores": midrank_z,
              "original scale": lambda x: (x - x.mean()) / x.std()}
    if var.startswith("Earthworms"):
        scales["log(1 + count)"] = lambda x: np.log1p(x)
        scales["square root of count"] = lambda x: np.sqrt(x)
    for sname, fn in scales.items():
        y = np.full(len(v), np.nan); y[ok] = fn(v[ok])
        r = dict(variable=var, scale=sname, **shares(y, ok))
        rg = pd.Series(cl[ok]); g = cl[ok]
        r["ICC1_collector"] = icc1(y[ok], g)
        rows.append(r)
        print(var, sname, {k: round(float(x), 3) for k, x in r.items() if isinstance(x, (float, np.floating))}, flush=True)
io.save_table(pd.DataFrame(rows), "T09_variance_share_by_scale.csv")
