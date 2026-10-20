import _bootstrap  # noqa
import warnings, json
import numpy as np, pandas as pd
from scipy import stats
import statsmodels.api as sm
from echo_ro import config as cfg, io
from echo_ro.stats_utils import grid_ids, blom_inverse_normal
warnings.filterwarnings("ignore")
d = io.load()
d["cl"] = pd.factorize(d.collected_by)[0]
d["cell25"] = grid_ids(d.x_km.values, d.y_km.values, 25.0)
d["block100"] = grid_ids(d.x_km.values, d.y_km.values, 100.0)
d["one"] = 1
V = {"pH": "soil_ph_field", "Organic matter": "soil_organic_matter", "Earthworms (counts)": "earthworms_count"}
COMPS = ["block100", "cell25", "collector"]
def fit(df, comps, fixed="y ~ 1", reml=True):
    vc = {k: f"0 + C({k})" for k in comps}
    m = sm.MixedLM.from_formula(fixed, groups="one", vc_formula=vc, re_formula="0", data=df)
    return m.fit(reml=reml, method=["powell"])
def comps_of(res): return dict(zip(list(res.model.exog_vc.names), np.atleast_1d(res.vcomp)))
def partition(sub, label, extra_fixed=None):
    rows = []
    for var, col in V.items():
        s = sub[~sub[col].isna()].copy()
        s["y"] = blom_inverse_normal(s[col].values.astype(float))
        s["collector"] = s["cl"]
        fixed = "y ~ 1" if not extra_fixed else f"y ~ {extra_fixed}"
        full = fit(s, COMPS, fixed)
        vcv = comps_of(full); tot = sum(vcv.values()) + full.scale
        row = dict(scenario=label, variable=var, n=len(s), collectors=s.collector.nunique(),
                   **{f"share_{k}": vcv[k] / tot for k in COMPS}, share_residual=full.scale / tot)
        red = fit(s, ["block100", "cell25"], fixed)
        lr = max(2 * (full.llf - red.llf), 0); row["LR_collector"] = lr; row["p_collector"] = 0.5 * stats.chi2.sf(lr, 1)
        rows.append(row); print(label, var, {k: round(float(x), 3) if isinstance(x, (float, np.floating)) else x for k, x in row.items()}, flush=True)
    return rows
out = []
out += partition(d, "all samples (reference)")
# A. suspected group accounts: >=5 samples on a single calendar day
per = d.groupby(["collected_by", "date"]).size().groupby("collected_by").max()
susp5 = per[per >= 5].index
susp8 = per[per >= 8].index
print("susp5", len(susp5), d.collected_by.isin(susp5).sum(), "susp8", len(susp8), d.collected_by.isin(susp8).sum())
out += partition(d[~d.collected_by.isin(susp5)], "excluding accounts with >=5 samples in one day")
out += partition(d[~d.collected_by.isin(susp8)], "excluding accounts with >=8 samples in one day")
# B. seasonality
d["season"] = pd.cut(pd.to_datetime(d.date).dt.month + 12 * (pd.to_datetime(d.date).dt.year - 2025), bins=[0, 12.5, 15.5, 16.5, 99], labels=["Nov-Dec 2025", "Jan-Mar 2026", "Apr 2026", "May-Oct 2026"])
print(d.season.value_counts())
ew = d.copy()
ew = ew[~ew.earthworms_count.isna()]
desc = ew.groupby("season").earthworms_count.agg(n="size", mean="mean", median="median", zero=lambda s: (s == 0).mean()).reset_index()
print(desc)
# season fixed effect in the crossed model for earthworms (and pH/OM for reference): LRT with ML
rows = []
for var, col in V.items():
    s = d[~d[col].isna()].copy(); s["y"] = blom_inverse_normal(s[col].values.astype(float)); s["collector"] = s["cl"]
    f1 = fit(s, COMPS, "y ~ C(season)", reml=False); f0 = fit(s, COMPS, "y ~ 1", reml=False)
    lr = max(2 * (f1.llf - f0.llf), 0); p = stats.chi2.sf(lr, 3)
    rows.append(dict(variable=var, LR_season=lr, df=3, p_season=p)); print("season FE", var, round(lr, 2), p, flush=True)
seas = pd.DataFrame(rows)
out += partition(d, "season (4 classes) as fixed effect", extra_fixed="C(season)")
out += partition(d[d.season == "Apr 2026"], "April 2026 only")
pd.DataFrame(out).to_csv("../sens_variance.csv", index=False)
desc.to_csv("../sens_earthworm_by_season.csv", index=False); seas.to_csv("../sens_season_FE.csv", index=False)
mm = d.groupby("month").earthworms_count.agg(n="size", mean="mean", zero=lambda s: (s == 0).mean()).reset_index(); mm.to_csv("../sens_earthworm_by_month.csv", index=False)
json.dump(dict(susp5=[len(susp5), int(d.collected_by.isin(susp5).sum())], susp8=[len(susp8), int(d.collected_by.isin(susp8).sum())]), open("../sens_flags.json", "w"))
