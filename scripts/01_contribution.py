"""Step 1 (O1) - contribution per collector, concentration and temporal pattern.

Methods: descriptive statistics; Wilson 95 % CI for the share of single-sample collectors;
Gini coefficient with percentile bootstrap CI (resampling collectors); Lorenz curve;
Hill numbers (N0, N1, N2) as effective numbers of contributors; share of samples
contributed by the most active collectors; number of distinct active days; spatial
footprint (maximum distance among a collector's samples; distances carry the +-2 km error).
"""
import _bootstrap  # noqa: F401
import numpy as np, pandas as pd
from statsmodels.stats.proportion import proportion_confint
from echo_ro import config as C, io
from echo_ro.stats_utils import gini, hill_numbers, haversine_matrix

rng = io.rng("m01")
d = io.load()
N = len(d)
c = d.groupby("collected_by").size().sort_values(ascending=False)
n = c.values
G = len(n)
R = {"n_samples": N, "n_collectors": G}
R["per_collector"] = dict(mean=n.mean(), median=float(np.median(n)), q25=float(np.quantile(n, .25)),
                          q75=float(np.quantile(n, .75)), max=int(n.max()), sd=n.std(ddof=1))
single = int((n == 1).sum())
lo, hi = proportion_confint(single, G, method="wilson")
R["single"] = dict(n=single, pct=single / G * 100, wilson_ci_pct=[lo * 100, hi * 100],
                   share_of_samples_pct=single / N * 100)
bs = [gini(rng.choice(n, G, replace=True)) for _ in range(C.N_BOOT)]
R["gini"] = dict(value=gini(n), ci=list(np.quantile(bs, [.025, .975])))
cs = np.cumsum(n) / N
for p in [1, 5, 10, 20]:
    k = int(np.ceil(G * p / 100))
    R[f"top{p}pct"] = dict(k=k, share_pct=n[:k].sum() / N * 100)
R["collectors_for_50pct_of_samples"] = int(np.searchsorted(cs, .5) + 1)
R["collectors_for_80pct_of_samples"] = int(np.searchsorted(cs, .8) + 1)
R["hill"] = hill_numbers(n)
R["hill"]["evenness_N2_over_N0"] = R["hill"]["N2"] / G

g = d.groupby("collected_by")
days = g["date"].nunique()
span = g["dt"].agg(lambda s: (s.max() - s.min()).days)
R["activity"] = dict(median_active_days=float(days.median()), multi_day_pct=(days > 1).mean() * 100,
                     multi_day_pct_among_ge2=(days[c >= 2] > 1).mean() * 100,
                     median_span_days_multiday=float(span[days > 1].median()),
                     span_gt30d_pct=(span > 30).mean() * 100,
                     max_samples_per_collector_day=int(d.groupby(["collected_by", "date"]).size().max()))
D = haversine_matrix(d.lat.values, d.lon.values)
fp = [D[np.ix_(ix, ix)].max() for ix in d.groupby("collected_by").indices.values() if len(ix) >= 2]
fp = np.array(fp)
R["footprint_ge2_samples"] = dict(n_collectors=len(fp), median_km=float(np.median(fp)),
                                  le5km_pct=(fp <= 5).mean() * 100, le20km_pct=(fp <= 20).mean() * 100,
                                  gt50km_pct=(fp > 50).mean() * 100)
top10 = set(c.index[:int(np.ceil(G * .10))])
apr = d.month == "2026-04"
R["campaign"] = dict(april_samples=int(apr.sum()), april_pct=apr.mean() * 100,
                     collectors_active_april=int(d[apr].collected_by.nunique()),
                     collectors_only_april=int(g.month.agg(lambda s: set(s) == {"2026-04"}).sum()),
                     top10pct_share_in_april=d[apr].collected_by.isin(top10).mean() * 100,
                     top10pct_share_other_months=d[~apr].collected_by.isin(top10).mean() * 100)
io.save_json(R, "T01_contribution_summary.json")
io.save_table(c.rename("n_samples").reset_index(), "T01b_samples_per_collector.csv")
m = d.groupby("month").agg(samples=("sample_id", "size"), collectors=("collected_by", "nunique")).reset_index()
io.save_table(m, "T01c_samples_per_month.csv")
print("gini", R["gini"], "single", R["single"]["pct"])
