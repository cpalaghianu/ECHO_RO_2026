"""Step 2 (O1) - completeness, domain checks and plausibility of the field indicators.

Methods
- Completeness / domain: missing values, duplicate identifiers, admissible value ranges.
- Prevalence of selected classes: Wilson 95 % CI (assumes independent samples) and percentile
  cluster bootstrap CI that resamples collectors (accounts for several samples per collector).
- Digit preference (heaping) of pH: observed frequency of whole-unit readings relative to the mean
  frequency of the two adjacent half-unit readings (ratio 1 = no preference); cluster-bootstrap CI.
  Heterogeneity of the whole/half-unit share among collectors with >= 10 samples: chi-square
  statistic with a Monte-Carlo permutation p-value.
- Earthworm counts: Poisson dispersion test (sum (y-mean)^2/mean ~ chi2_{n-1}); intercept-only
  negative-binomial (NB2) fit and expected number of zeros under Poisson and NB;
  heaping at 5 and 10 relative to the geometric interpolation of the neighbouring counts
  (Poisson tail probability, exploratory).
- Concentration of extreme classes (pH <= 3.5, pH >= 9, pollutants >= 12) by collector.
- Plausibility of pH <= 3.5: median pH of the five nearest eligible neighbours (>= 2 km, other
  collector) compared with that of samples with pH 4.0-4.5 (Mann-Whitney U).
"""
import _bootstrap  # noqa: F401
import numpy as np, pandas as pd
from scipy import stats
from statsmodels.stats.proportion import proportion_confint
import statsmodels.api as sm
from echo_ro import config as C, io
from echo_ro.stats_utils import haversine_matrix

rng = io.rng("m02")
d = io.load()
N = len(d)
R = {}
cols = ["sample_id", "sampling_datetime_utc", "lat", "lon", "soil_ph_field", "soil_organic_matter",
        "earthworms_count", "pollutants_count", "soil_structure_en", "soil_texture_en", "collected_by"]
R["missing"] = {c: int(d[c].isna().sum()) for c in cols}
R["duplicated_ids"] = int(d.sample_id.duplicated().sum())
R["ph_on_half_unit_grid"] = bool(((d.soil_ph_field * 2) % 1 == 0).all())
R["ph_range"] = [d.soil_ph_field.min(), d.soil_ph_field.max()]
R["counts_are_nonnegative_integers"] = bool(((d.earthworms_count >= 0) & (d.pollutants_count >= 0)).all())

# --- prevalence with Wilson and cluster-bootstrap CIs
props = {
    "pH <= 4.5": d.soil_ph_field <= 4.5, "pH <= 3.5": d.soil_ph_field <= 3.5, "pH >= 7.5": d.soil_ph_field >= 7.5,
    "OM lowest two classes (1.5, 2)": d.soil_organic_matter.isin([1.5, 2]),
    "Structure friable": d.soil_structure_en == "friable",
    "Structure compact or high compact": d.soil_structure_en.isin(["compact", "high compact"]),
    "Earthworms = 0": d.earthworms_count == 0, "Pollutants >= 1": d.pollutants_count >= 1,
    "Texture = other": d.soil_texture_en == "other"}
rows = []
for k, v in props.items():
    ok = d.soil_organic_matter.notna() if k.startswith("OM") else pd.Series(True, index=d.index)
    sub = d[ok]
    ind = v[ok].values.astype(float)
    u, inv = np.unique(sub.collected_by.values, return_inverse=True)
    num, den = np.bincount(inv, weights=ind), np.bincount(inv).astype(float)
    bs = []
    for _ in range(C.N_BOOT):
        s = rng.integers(0, len(u), len(u))
        bs.append(num[s].sum() / den[s].sum())
    w = proportion_confint(int(ind.sum()), int(ok.sum()), method="wilson")
    rows.append(dict(indicator=k, x=int(ind.sum()), n=int(ok.sum()), pct=ind.mean() * 100,
                     wilson_lo=w[0] * 100, wilson_hi=w[1] * 100,
                     boot_lo=np.quantile(bs, .025) * 100, boot_hi=np.quantile(bs, .975) * 100))
io.save_table(pd.DataFrame(rows), "T02_prevalence_ci.csv")

# --- pH heaping
u, inv = np.unique(d.collected_by.values, return_inverse=True)
G = len(u)
ix_by = [np.where(inv == g)[0] for g in range(G)]
phv = d.soil_ph_field.values
vc = d.soil_ph_field.value_counts().sort_index()
whole = [v for v in vc.index if v == int(v)]
half = [v for v in vc.index if v != int(v)]
inner = [w for w in whole if (w - .5 in vc.index) and (w + .5 in vc.index)]
io.save_table(vc.rename("n").rename_axis("ph").reset_index(), "T02b_ph_frequency.csv")

def heap_ratio(values):
    c = pd.Series(values).value_counts()
    o = sum(c.get(w, 0) for w in inner)
    e = sum((c.get(w - .5, 0) + c.get(w + .5, 0)) / 2 for w in inner)
    return o / e
bs = [heap_ratio(phv[np.concatenate([ix_by[g] for g in rng.integers(0, G, G)])]) for _ in range(C.N_BOOT // 2)]
R["ph_whole_share_pct"] = vc[whole].sum() / N * 100
R["ph_whole_share_expected_uniform_pct"] = len(whole) / len(vc) * 100
R["ph_heaping_ratio"] = dict(value=heap_ratio(phv), ci=list(np.quantile(bs, [.025, .975])))
big = [g for g in range(G) if len(ix_by[g]) >= 10]
tab = np.array([[np.isin(phv[ix_by[g]], whole).sum(), np.isin(phv[ix_by[g]], half).sum()] for g in big])
sizes = tab.sum(1)
lab = np.concatenate([np.r_[np.ones(a), np.zeros(b)] for a, b in tab])
def chi_stat(l):
    pos, o = 0, []
    for s in sizes:
        o.append(l[pos:pos + s].sum()); pos += s
    o = np.array(o); e = sizes * lab.mean()
    return (((o - e) ** 2 / e) + (((sizes - o) - (sizes - e)) ** 2 / (sizes - e))).sum()
obs_c = chi_stat(lab)
cnt = sum(chi_stat(rng.permutation(lab)) >= obs_c for _ in range(C.N_PERM))
R["ph_heaping_heterogeneity"] = dict(n_collectors=len(big), n_samples=int(sizes.sum()), chi2=obs_c,
                                     p_monte_carlo=(cnt + 1) / (C.N_PERM + 1),
                                     half_unit_share_range_pct=[float((tab[:, 1] / sizes).min() * 100),
                                                                float((tab[:, 1] / sizes).max() * 100)])

# --- extreme classes
def conc(mask):
    byc = d[mask].groupby("collected_by").size().sort_values(ascending=False)
    return dict(n=int(mask.sum()), collectors=int(len(byc)), top_collector_n=int(byc.iloc[0]),
                top3_share_pct=float(byc.iloc[:3].sum() / mask.sum() * 100), months=d[mask].month.value_counts().to_dict())
R["ph_le_3.5"] = conc(d.soil_ph_field <= 3.5)
R["ph_ge_9"] = conc(d.soil_ph_field >= 9)
R["texture_other"] = conc(d.soil_texture_en == "other")
R["pollutants_ge_12"] = d[d.pollutants_count >= 12][["sample_id", "collected_by", "pollutants_count", "soil_ph_field"]].to_dict("records")

# --- earthworms
y = d.earthworms_count.values
mu, var = y.mean(), y.var(ddof=1)
chi = ((y - mu) ** 2 / mu).sum()
nb = sm.NegativeBinomial(y, np.ones(N)).fit(disp=0, maxiter=200)
a, m_ = float(nb.params[1]), float(np.exp(nb.params[0]))
W = dict(mean=mu, variance=var, variance_to_mean=var / mu, zeros=int((y == 0).sum()), zero_pct=(y == 0).mean() * 100,
         max=int(y.max()), poisson_dispersion_chi2=chi, poisson_dispersion_p=float(stats.chi2.sf(chi, N - 1)),
         nb_alpha=a, nb_mean=m_, zeros_expected_poisson=float(N * np.exp(-mu)),
         zeros_expected_nb=float(N * (1 + a * m_) ** (-1 / a)))
c = pd.Series(y).value_counts()
for k, (lo, hi) in {5: (4, 6), 10: (9, 11)}.items():
    e = np.sqrt(c.get(lo, 0) * c.get(hi, 0)); o = c.get(k, 0)
    W[f"heaping_{k}"] = dict(observed=int(o), expected_geometric=float(e), ratio=float(o / e),
                             p_poisson_upper=float(stats.poisson.sf(o - 1, e)))
R["earthworms"] = W
kk = np.arange(0, int(y.max()) + 1)
from scipy.stats import nbinom
r_ = 1.0 / a
obs_f = np.array([(y == v).sum() for v in kk])
exp_nb = N * nbinom.pmf(kk, r_, r_ / (r_ + m_))
exp_po = N * stats.poisson.pmf(kk, mu)
io.save_table(pd.DataFrame(dict(count=kk, observed=obs_f, expected_nb=exp_nb, expected_poisson=exp_po)), "T02c_earthworm_fit.csv")

p = d.pollutants_count
R["pollutants"] = dict(zero_pct=(p == 0).mean() * 100, ge1=int((p >= 1).sum()), ge12=int((p >= 12).sum()))

# --- plausibility of very low pH against neighbours (eligible = >= 2 km and different collector)
D = haversine_matrix(d.lat.values, d.lon.values)
np.fill_diagonal(D, np.inf)
D[D < C.MIN_PAIR_DISTANCE_KM] = np.inf
cl = pd.factorize(d.collected_by)[0]
D[cl[:, None] == cl[None, :]] = np.inf
nbi = np.argsort(D, 1)[:, :5]
t = pd.DataFrame({"ph": d.soil_ph_field.values, "nb_median": np.median(phv[nbi], 1)})
g = t.groupby("ph").nb_median.agg(n="count", median="median", q25=lambda s: s.quantile(.25), q75=lambda s: s.quantile(.75))
io.save_table(g.reset_index(), "T02d_neighbour_median_by_ph.csv")
lo_, mid_ = t[t.ph <= 3.5].nb_median, t[(t.ph >= 4) & (t.ph <= 4.5)].nb_median
mw = stats.mannwhitneyu(lo_, mid_, alternative="two-sided")
R["ph_le_3.5_neighbour_check"] = dict(median_nb_ph_le_3_5=float(lo_.median()), median_nb_ph_4_to_4_5=float(mid_.median()),
                                      mannwhitney_U=float(mw.statistic), p=float(mw.pvalue),
                                      share_nb_ge_5_5_ph_le_3_5=float((lo_ >= 5.5).mean()),
                                      share_nb_ge_5_5_ph_4_to_4_5=float((mid_ >= 5.5).mean()))
io.save_json(R, "T02_quality_summary.json")
print(R["ph_heaping_ratio"], R["earthworms"]["heaping_5"], R["ph_le_3.5_neighbour_check"])
