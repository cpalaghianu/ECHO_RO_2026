"""Step 6 - publication figures (matplotlib only), generated from the tables in outputs/tables.

Fig. 1  Network: Lorenz curve, contributions per collector, seasonal pattern of sampling.
Fig. 2  Data quality: pH digit preference, earthworm counts vs. fitted distributions, prevalence of classes.
Fig. 3  Spatial autocorrelation correlograms, nominal-variable agreement and variance partition.
Fig. 4  Neighbour consistency (naive vs. corrected) and collector effect.
Output: outputs/figures/FigN_*.pdf (vector, embedded TrueType) and .png (300 dpi).
"""
import _bootstrap  # noqa: F401
import json
import numpy as np, pandas as pd
import matplotlib.pyplot as plt
from matplotlib.patches import Patch
from matplotlib.lines import Line2D
from echo_ro import config as cfg, style as S

S.apply()
T = lambda f: pd.read_csv(cfg.TABLES / f)
J = lambda f: json.load(open(cfg.TABLES / f))
BANDS = ["2-5 km", "5-10 km", "10-20 km", "20-50 km", "50-100 km", "100-200 km"]
BL = ["2–5", "5–10", "10–20", "20–50", "50–100", "100–200"]

# ---------------------------------------------------------------- Fig 1
c = T("T01b_samples_per_collector.csv").n_samples.values
G, N = len(c), c.sum()
mo = T("T01c_samples_per_month.csv")
R1 = J("T01_contribution_summary.json")
fig, ax = plt.subplots(1, 3, figsize=(S.W2, 55 * S.MM), gridspec_kw=dict(wspace=0.45))
xs = np.r_[0, np.arange(1, G + 1) / G]; ys = np.r_[0, np.cumsum(np.sort(c)) / N]
ax[0].plot([0, 1], [0, 1], color=S.GREY, lw=0.8, ls=":")
ax[0].plot(xs, ys, color=S.BLUE, lw=1.4)
ax[0].fill_between(xs, xs, ys, color=S.BLUE, alpha=0.10, lw=0)
ax[0].set_xlabel("Cumulative share of collectors"); ax[0].set_ylabel("Cumulative share of samples")
ax[0].text(0.05, 0.9, f"Gini = {R1['gini']['value']:.2f}", transform=ax[0].transAxes)
S.panel_label(ax[0], "a")
bins = [(1, 1), (2, 2), (3, 3), (4, 4), (5, 5), (6, 10), (11, 20), (21, 99)]
cnt = [int(((c >= lo) & (c <= hi)).sum()) for lo, hi in bins]
lbl = [str(lo) if lo == hi else (f"{lo}–{hi}" if hi < 99 else f">{lo-1}") for lo, hi in bins]
ax[1].bar(range(len(bins)), cnt, color=S.BLUE, width=0.75)
ax[1].set_yscale("log"); ax[1].set_ylim(1, 400); ax[1].set_xticks(range(len(bins))); ax[1].set_xticklabels(lbl, rotation=45, ha="right")
ax[1].set_xlabel("Samples per collector"); ax[1].set_ylabel("Number of collectors (log scale)")
ax[1].text(0.97, 0.95, f"{R1['single']['pct']:.0f}% of collectors\ncontributed one sample", transform=ax[1].transAxes, ha="right", va="top")
S.panel_label(ax[1], "b")
x = np.arange(len(mo)); lab = [pd.Timestamp(m + "-01").strftime("%b %y") for m in mo.month]
ax[2].bar(x, mo.samples, color=S.LIGHT, width=0.75, label="Samples")
ax[2].plot(x, mo.collectors, color=S.BLUE, marker="o", lw=1.2, label="Active collectors")
ax[2].set_xticks(x); ax[2].set_xticklabels(lab, rotation=60, ha="right"); ax[2].set_ylabel("Count")
ax[2].legend(loc="upper right", handlelength=1.4); S.panel_label(ax[2], "c")
S.save(fig, "Fig1_network", cfg.FIGURES)

# ---------------------------------------------------------------- Fig 2
ph = T("T02b_ph_frequency.csv"); ew = T("T02c_earthworm_fit.csv"); pr = T("T02_prevalence_ci.csv")
fig = plt.figure(figsize=(S.W2, 100 * S.MM))
gs = fig.add_gridspec(2, 2, hspace=0.95, wspace=0.35, height_ratios=[1, 1.15])
ax = [fig.add_subplot(gs[0, 0]), fig.add_subplot(gs[0, 1]), fig.add_subplot(gs[1, :])]
whole = (ph.ph % 1 == 0)
ax[0].bar(ph.ph[whole], ph.n[whole], width=0.38, color=S.BLUE, label="Whole units")
ax[0].bar(ph.ph[~whole], ph.n[~whole], width=0.38, color="white", edgecolor=S.GREY, hatch="////", lw=0.6, label="Half units")
ax[0].set_xlabel("Field pH"); ax[0].set_ylabel("Number of samples")
ax[0].legend(loc="lower center", bbox_to_anchor=(0.5, 1.0), ncol=2, handlelength=1.4, columnspacing=1.0)
ax[0].set_xticks(range(3, 10)); S.panel_label(ax[0], "a")
k = ew[ew["count"] <= 20]
ax[1].bar(k["count"], k.observed, color=S.LIGHT, width=0.8, label="Observed")
ax[1].plot(k["count"], k.expected_nb, color=S.OCHRE, marker="s", ls="--", mec=S.INK, mew=0.5, ms=3.2, label="Negative binomial")
ax[1].plot(k["count"], k.expected_poisson, color=S.TEAL, marker="^", ls=":", mec=S.INK, mew=0.5, ms=3.2, label="Poisson")
ax[1].set_xlabel("Earthworms per sample"); ax[1].set_ylabel("Number of samples")
ax[1].legend(loc="center right", handlelength=1.4)
ax[1].set_xticks(range(0, 21, 5)); S.panel_label(ax[1], "b")
pr = pr.iloc[::-1].reset_index(drop=True); y = np.arange(len(pr))
ax[2].grid(axis="x"); ax[2].grid(axis="y", visible=False)
ax[2].hlines(y + 0.14, pr.wilson_lo, pr.wilson_hi, color=S.GREY, lw=1.0)
ax[2].plot(pr.pct, y + 0.14, "s", mfc="white", mec=S.GREY, ms=3.5)
ax[2].hlines(y - 0.14, pr.boot_lo, pr.boot_hi, color=S.BLUE, lw=1.4)
ax[2].plot(pr.pct, y - 0.14, "o", color=S.BLUE, ms=3.5)
ax[2].set_yticks(y); ax[2].set_yticklabels([t.replace("<=", "≤").replace(">=", "≥") for t in pr.indicator])
ax[2].set_xlabel("Share of samples (%) with 95% confidence interval")
ax[2].legend(handles=[Line2D([], [], color=S.GREY, marker="s", mfc="white", lw=1, label="Wilson (independent samples)"),
                      Line2D([], [], color=S.BLUE, marker="o", lw=1.4, label="Cluster bootstrap (collectors resampled)")],
             loc="upper right", handlelength=1.4)
S.panel_label(ax[2], "c")
S.save(fig, "Fig2_data_quality", cfg.FIGURES)

# ---------------------------------------------------------------- Fig 3
sp = T("T03_spatial_autocorrelation.csv"); vp = T("T05b_variance_partition_crossed.csv")
fig, ax = plt.subplots(2, 3, figsize=(S.W2, 108 * S.MM), gridspec_kw=dict(hspace=0.8, wspace=0.55))
xb = np.arange(len(BANDS))
def corr(a, var, ylabel, col="statistic", nullname="block100", letter=""):
    g = sp[(sp.variable == var) & (sp["null"] == "global")]
    b = sp[(sp.variable == var) & (sp["null"] == nullname)]
    nb = b[b.pairs == "all"].set_index("band").reindex(BANDS)
    ok = nb.null_lo.notna().values
    lo_, hi_ = nb.null_lo, nb.null_hi
    if col == "kappa_s":       # express the permutation band on the kappa scale
        lo_ = (nb.null_lo - nb.null_mean) / (1 - nb.null_mean); hi_ = (nb.null_hi - nb.null_mean) / (1 - nb.null_mean)
    a.fill_between(xb[ok], lo_[ok], hi_[ok], color=S.LIGHT, alpha=0.7, lw=0, hatch="////", edgecolor="white")
    for (pr_, st) in [("all", S.SERIES[0]), ("between", S.SERIES[1])]:
        s = g[g.pairs == pr_].set_index("band").reindex(BANDS)
        sb = b[b.pairs == pr_].set_index("band").reindex(BANDS)
        a.plot(xb, s[col], color=st["color"], ls=st["ls"], lw=1.2)
        sig = (sb.q_BH < 0.05).reindex(BANDS).fillna(False).values
        for i in range(len(BANDS)):
            filled = bool(sig[i])
            a.plot(xb[i], s[col].iloc[i], st["marker"], ms=4.5, mec=S.INK if st["color"] == S.OCHRE else st["color"], mew=0.8,
                   mfc=st["color"] if filled else "white")
    a.axhline(0, color=S.INK, lw=0.6)
    a.set_xticks(xb); a.set_xticklabels(BL, rotation=45, ha="right"); a.set_xlabel("Distance band (km)"); a.set_ylabel(ylabel)
    S.panel_label(a, letter)
corr(ax[0, 0], "pH", "Moran's I", letter="a"); ax[0, 0].set_title("pH", loc="left", pad=3)
corr(ax[0, 1], "Organic matter", "Moran's I", letter="b"); ax[0, 1].set_title("Organic matter", loc="left", pad=3)
corr(ax[0, 2], "Earthworms (counts)", "Moran's I", letter="c"); ax[0, 2].set_title("Earthworm counts", loc="left", pad=3)
corr(ax[1, 0], "Structure", "Spatial kappa", col="kappa_s", letter="d"); ax[1, 0].set_title("Structure (nominal)", loc="left", pad=3)
corr(ax[1, 1], "Texture (excl. other)", "Spatial kappa", col="kappa_s", letter="e"); ax[1, 1].set_title("Texture (nominal)", loc="left", pad=3)
leg = [Line2D([], [], color=S.BLUE, marker="o", ls="-", lw=1.2, label="All pairs"),
       Line2D([], [], color=S.OCHRE, marker="s", ls="--", mec=S.INK, lw=1.2, label="Different collectors"),
       Patch(facecolor=S.LIGHT, hatch="////", edgecolor="white", label="95% null, 100-km blocks"),
       Line2D([], [], color=S.INK, marker="o", mfc="white", ls="", label="Filled: q < 0.05 vs. block null")]
ax[0, 0].legend(handles=leg, loc="lower center", bbox_to_anchor=(1.75, 1.18), ncol=4, columnspacing=1.2, handlelength=1.6)
a = ax[1, 2]
vv = vp.set_index("variable").loc[["Earthworms (counts)", "Organic matter", "pH"]]
place = vv.share_block100 + vv.share_cell25
yy = np.arange(3)
a.barh(yy, vv.share_collector * 100, color=S.BLUE, height=0.55, label="Collector")
a.barh(yy, place * 100, left=vv.share_collector * 100, color=S.TEAL, height=0.55, label="Place (25-km cell + 100-km block)")
a.barh(yy, vv.share_residual * 100, left=(vv.share_collector + place) * 100, color="white", edgecolor=S.GREY, hatch="////", lw=0.6, height=0.55, label="Residual")
a.set_yticks(yy); a.set_yticklabels(["Earthworms", "Org. matter", "pH"]); a.set_xlim(0, 100)
a.set_xlabel("Share of variance (%)"); a.grid(axis="x"); a.grid(axis="y", visible=False)
a.legend(loc="upper center", bbox_to_anchor=(0.5, -0.36), ncol=1, handlelength=1.4); S.panel_label(a, "f")
S.save(fig, "Fig3_spatial_structure", cfg.FIGURES)

# ---------------------------------------------------------------- Fig 4
kn = T("T04_knn_consistency.csv"); op = T("T04b_collector_effect.csv")
fig, ax = plt.subplots(1, 2, figsize=(S.W2, 52 * S.MM), gridspec_kw=dict(width_ratios=[1, 1], wspace=0.75))
vars_ = ["pH", "Organic matter", "Earthworms (counts)"]; yy = np.arange(len(vars_))[::-1]
a = ax[0]; a.grid(axis="x"); a.grid(axis="y", visible=False)
for mode, d0, off, st in [("naive", 0.0, 0.16, S.SERIES[1]), ("corrected", 2.0, -0.16, S.SERIES[0])]:
    s = kn[(kn["mode"] == mode) & (kn.d0_km == d0) & (kn.k == 5)].set_index("variable").loc[vars_]
    a.hlines(yy + off, s.ci_blockboot_lo, s.ci_blockboot_hi, color=st["color"], lw=1.4)
    a.plot(s.rho, yy + off, st["marker"], color=st["color"], mec=S.INK if st["color"] == S.OCHRE else st["color"], ms=4.5,
           label="All other samples" if mode == "naive" else "Eligible neighbours (≥ 2 km, other collector)")
a.axvline(0, color=S.INK, lw=0.6); a.set_yticks(yy); a.set_yticklabels(["pH", "Organic matter", "Earthworms"])
a.set_xlabel("Spearman ρ with median of 5 neighbours (95% block-bootstrap CI)"); a.legend(loc="upper center", bbox_to_anchor=(0.45, -0.28), handlelength=1.4)
S.panel_label(a, "a")
a = ax[1]; a.grid(axis="x"); a.grid(axis="y", visible=False)
op = op.assign(ratio=op.mean_dissimilarity_same / op.mean_dissimilarity_other).sort_values("ratio", ascending=False)
yy = np.arange(len(op))[::-1]
a.hlines(yy, op.ratio, 1, color=S.LIGHT, lw=1.2)
a.scatter(op.ratio, yy, s=22, marker="o", facecolor=[S.BLUE if p < 0.05 else "white" for p in op.p_perm], edgecolor=S.BLUE, zorder=3)
a.axvline(1, color=S.INK, lw=0.6); a.set_yticks(yy); a.set_yticklabels(op.variable)
a.set_xlabel("Dissimilarity ratio, same / different collector"); a.set_xlim(0, 1.15)
a.text(0.02, -0.38, "Filled: permutation p < 0.05", transform=a.transAxes, color=S.INK)
S.panel_label(a, "b")
S.save(fig, "Fig4_consistency", cfg.FIGURES)
print("figures written to", cfg.FIGURES)
