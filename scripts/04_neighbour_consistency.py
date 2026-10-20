"""Step 4 (O4) - plausible consistency of neighbouring samples and collector (operator) effect.

A. k-nearest-neighbour (kNN) consistency
   For each sample the median value of its k nearest neighbours is computed; consistency is the
   Spearman rank correlation between the sample value and its neighbourhood median.
   'naive'     : all other samples are eligible neighbours;
   'corrected' : neighbours must be >= d0 km away (d0 = 2 or 5 km; coordinate error) AND
                 collected by a different person (removes operator and same-site dependence).
   Uncertainty: percentile CIs from (i) a cluster bootstrap resampling collectors and (ii) a spatial
   block bootstrap resampling 100-km blocks (the latter also accounts for spatial dependence and
   is the more conservative interval).
   Permutation tests (k = 5, d0 = 2 km, neighbour sets fixed): global shuffling and shuffling
   within 100-km blocks (the latter tests whether local coherence exceeds the regional gradient).
B. Collector effect conditional on location
   Among sample pairs 2-25 km apart located in the same 25-km grid cell, the mean dissimilarity of
   pairs from the same collector is compared with pairs from different collectors
   (T = mean_same - mean_different; dissimilarity = |rank difference| / n for ordinal variables,
   0/1 mismatch for nominal ones). p-value: collector labels are permuted within 25-km cells
   (keeps spatial composition, breaks collector identity).
Usage: python scripts/04_neighbour_consistency.py [--om complete_case|imputed]
"""
import _bootstrap  # noqa: F401
import argparse
import numpy as np, pandas as pd
from scipy import stats
from echo_ro import config as C, io
from echo_ro.stats_utils import haversine_matrix, grid_ids, groups_from_ids

ap = argparse.ArgumentParser()
ap.add_argument("--om", default="complete_case", choices=["complete_case", "imputed"])
args = ap.parse_args()
suffix = "" if args.om == "complete_case" else "_omimputed"
rng = io.rng("m04")
d = io.load(om=args.om)
N = len(d)
D = haversine_matrix(d.lat.values, d.lon.values)
cl = pd.factorize(d.collected_by)[0]
V = {"pH": d.soil_ph_field.values.astype(float), "Organic matter": d.soil_organic_matter.values.astype(float),
     "Earthworms (counts)": d.earthworms_count.values.astype(float)}
sp = lambda a, b: stats.spearmanr(a, b)[0]

rows = []
for var, v in V.items():
    idx = np.where(~np.isnan(v))[0]; m = len(idx); vm = v[idx]; clm = cl[idx]
    Dm = D[np.ix_(idx, idx)]
    blocks100 = groups_from_ids(grid_ids(d.x_km.values[idx], d.y_km.values[idx], 100.0))
    for mode, d0 in [("naive", 0.0)] + [("corrected", x) for x in C.KNN_D0]:
        for k in C.KNN_K:
            Dd = Dm.copy(); np.fill_diagonal(Dd, np.inf)
            if mode == "corrected":
                Dd[Dd < d0] = np.inf
                Dd[clm[:, None] == clm[None, :]] = np.inf
            nb = np.argsort(Dd, 1)[:, :k]; dist = np.take_along_axis(Dd, nb, 1)
            ok = np.isfinite(dist).all(1)
            own, nbm = vm[ok], np.median(vm[nb[ok]], 1)
            rho = sp(own, nbm)
            cc = clm[ok]; uc = np.unique(cc); ixs = [np.where(cc == c)[0] for c in uc]
            bs = []
            for _ in range(C.N_BOOT_KNN):
                s = rng.integers(0, len(uc), len(uc)); ii = np.concatenate([ixs[g] for g in s])
                bs.append(sp(own[ii], nbm[ii]))
            # spatial block bootstrap: 100-km blocks resampled (keeps collector and spatial dependence inside blocks)
            bk = grid_ids(d.x_km.values[idx], d.y_km.values[idx], 100.0)[ok]
            bix = [np.where(bk == b)[0] for b in np.unique(bk)]
            bsb = []
            for _ in range(C.N_BOOT_KNN):
                s_ = rng.integers(0, len(bix), len(bix)); ii = np.concatenate([bix[g] for g in s_])
                bsb.append(sp(own[ii], nbm[ii]))
            p_glob = p_blk = np.nan
            if mode == "corrected" and k == 5 and d0 == C.MIN_PAIR_DISTANCE_KM:
                okidx, nbo = np.where(ok)[0], nb[ok]
                def perm_p(blocks):
                    cnt = 0
                    for _ in range(C.N_PERM_KNN):
                        if blocks is None:
                            pv = vm[rng.permutation(m)]
                        else:
                            pi = np.arange(m)
                            for ib in blocks:
                                pi[ib] = ib[rng.permutation(len(ib))]
                            pv = vm[pi]
                        cnt += abs(sp(pv[okidx], np.median(pv[nbo], 1))) >= abs(rho)
                    return (cnt + 1) / (C.N_PERM_KNN + 1)
                p_glob, p_blk = perm_p(None), perm_p(blocks100)
            rows.append(dict(variable=var, mode=mode, d0_km=d0, k=k, n=int(ok.sum()), rho=rho,
                             ci_lo=np.quantile(bs, .025), ci_hi=np.quantile(bs, .975),
                             ci_blockboot_lo=np.quantile(bsb, .025), ci_blockboot_hi=np.quantile(bsb, .975),
                             median_mean_neighbour_distance_km=float(np.median(dist[ok].mean(1))),
                             p_perm_global=p_glob, p_perm_block100=p_blk))
    print(var, "kNN done", flush=True)
io.save_table(pd.DataFrame(rows), f"T04_knn_consistency{suffix}.csv")

# ---------- B. collector effect conditional on 25-km cells
cell25 = grid_ids(d.x_km.values, d.y_km.values, 25.0)
iu, ju = np.triu_indices(N, 1)
dd = D[iu, ju]
sel = (dd >= C.MIN_PAIR_DISTANCE_KM) & (dd < 25) & (cell25[iu] == cell25[ju])
iu, ju = iu[sel], ju[sel]
cb = groups_from_ids(cell25)
tests = {"pH": ("ord", d.soil_ph_field.values), "Organic matter": ("ord", d.soil_organic_matter.values),
         "Earthworms (counts)": ("ord", d.earthworms_count.values),
         "Structure": ("cat", pd.factorize(d.soil_structure_en)[0]),
         "Texture (excl. other)": ("cat", pd.factorize(d.soil_texture_en)[0]),
         "Pollutants present": ("cat", (d.pollutants_count.values > 0).astype(int)),
         "Earthworms present": ("cat", (d.earthworms_count.values > 0).astype(int))}
op = []
for name, (kind, v) in tests.items():
    ok = ~pd.isna(v) if kind == "ord" else np.ones(N, bool)
    if name.startswith("Texture"):
        ok = d.soil_texture_en.values != "other"
    use = ok[iu] & ok[ju]; a, b = iu[use], ju[use]
    if kind == "ord":
        rk = np.full(N, np.nan); rk[ok] = stats.rankdata(np.asarray(v, float)[ok])
        diff = np.abs(rk[a] - rk[b]) / ok.sum()
    else:
        diff = (v[a] != v[b]).astype(float)
    def T(lab):
        s = lab[a] == lab[b]
        return diff[s].mean() - diff[~s].mean(), int(s.sum())
    obs, ns = T(cl)
    nl = np.empty(C.N_PERM)
    for i in range(C.N_PERM):
        lab = cl.copy()
        for ib in cb:
            if len(ib) > 1:
                lab[ib] = cl[ib][rng.permutation(len(ib))]
        nl[i] = T(lab)[0]
    s0 = cl[a] == cl[b]
    op.append(dict(variable=name, pairs=len(a), same_collector_pairs=ns, mean_dissimilarity_same=diff[s0].mean(),
                   mean_dissimilarity_other=diff[~s0].mean(), T=obs,
                   p_perm=(1 + (np.abs(nl) >= abs(obs)).sum()) / (len(nl) + 1)))
    print(name, "operator test done", flush=True)
io.save_table(pd.DataFrame(op), f"T04b_collector_effect{suffix}.csv")
