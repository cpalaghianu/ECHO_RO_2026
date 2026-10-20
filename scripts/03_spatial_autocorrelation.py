"""Step 3 (O3, O4) - spatial autocorrelation of the field indicators by distance band.

Statistics
- Ordinal / count variables (pH, organic matter, earthworm counts): Moran's I on mid-ranks,
  with binary distance-band weights, row-standardised (Cliff & Ord 1981; Anselin 1995).
- Nominal variables (structure, texture excluding 'other', presence of earthworms, presence of
  pollutants): spatial agreement = mean share of neighbours (in the band) with the same class;
  spatial kappa = (observed - null mean) / (1 - null mean) (chance-corrected, cf. Cohen 1960).
Distance bands: 2-5, 5-10, 10-20, 20-50, 50-100, 100-200 km. Pairs closer than 2 km are never used,
because the provider displaced each coordinate by a random offset of up to 1 km (distance error
up to 2 km). Two pair sets: 'all' and 'between' (pairs of samples from different collectors only).
Inference (Monte-Carlo permutation, two-sided, (1 + #extreme) / (1 + B), B = N_PERM)
- 'global' null: values permuted over the whole country (tests any spatial structure,
  including the regional gradient);
- 'block50' / 'block100' nulls: values permuted only within 50- or 100-km grid blocks, which
  preserves the regional gradient and tests only the structure at scales smaller than the block
  (restricted permutation; cf. Besag & Diggle 1977; Legendre 1993). Applied to bands < 50 km (< block size).
Benjamini-Hochberg false-discovery-rate correction within each null model.
Usage: python scripts/03_spatial_autocorrelation.py [--om complete_case|imputed]
"""
import _bootstrap  # noqa: F401
import argparse, time
import numpy as np, pandas as pd
from scipy import stats, sparse
from statsmodels.stats.multitest import multipletests
from echo_ro import config as C, io
from echo_ro.stats_utils import haversine_matrix, grid_ids, groups_from_ids

ap = argparse.ArgumentParser()
ap.add_argument("--om", default="complete_case", choices=["complete_case", "imputed"])
args = ap.parse_args()
suffix = "" if args.om == "complete_case" else "_omimputed"
rng = io.rng("m03")
d = io.load(om=args.om)
N = len(d)
D = haversine_matrix(d.lat.values, d.lon.values)
cl = pd.factorize(d.collected_by)[0]
EDGES = C.DIST_EDGES_KM
NPERM, CH = C.N_PERM, 100
V = {
    "pH": ("ord", d.soil_ph_field.values.astype(float), np.ones(N, bool)),
    "Organic matter": ("ord", d.soil_organic_matter.values.astype(float), ~d.soil_organic_matter.isna().values),
    "Earthworms (counts)": ("ord", d.earthworms_count.values.astype(float), np.ones(N, bool)),
    "Earthworms present": ("cat", (d.earthworms_count.values > 0).astype(int), np.ones(N, bool)),
    "Pollutants present": ("cat", (d.pollutants_count.values > 0).astype(int), np.ones(N, bool)),
    "Structure": ("cat", pd.factorize(d.soil_structure_en)[0], np.ones(N, bool)),
    "Texture (excl. other)": ("cat", pd.factorize(d.soil_texture_en)[0], (d.soil_texture_en.values != "other")),
}

def make_perm(m, blocks, chunk):
    if blocks is None:
        return rng.random((chunk, m)).argsort(1)
    perm = np.empty((chunk, m), int)
    for ib in blocks:
        perm[:, ib] = ib[rng.random((chunk, len(ib))).argsort(1)]
    return perm

rows, t0 = [], time.time()
for name, (kind, v, valid) in V.items():
    idx = np.where(valid)[0]; m = len(idx); vm = v[idx]; clm = cl[idx]
    Dm = D[np.ix_(idx, idx)]
    iu, ju = np.triu_indices(m, 1)
    dd = Dm[iu, ju]; same = clm[iu] == clm[ju]
    if kind == "ord":
        z = stats.rankdata(vm); z = z - z.mean(); zz = (z ** 2).sum()
    else:
        z = pd.factorize(vm)[0]
    sets = {}
    for b in range(len(EDGES) - 1):
        sel = (dd >= EDGES[b]) & (dd < EDGES[b + 1])
        for mk in ["all", "between"]:
            s = sel if mk == "all" else sel & (~same)
            I = np.r_[iu[s], ju[s]]; J = np.r_[ju[s], iu[s]]
            if len(I) == 0:
                continue
            deg = np.bincount(I, minlength=m); w = 1.0 / deg[I]          # row-standardisation
            sets[(b, mk)] = dict(I=I, J=J, w=w, nb=int((deg > 0).sum()), npairs=int(s.sum()),
                                 W=sparse.csr_matrix((w, (I, J)), shape=(m, m)))
    def stat(Zp, st):
        if kind == "ord":      # Moran's I = (n_eff / S0) * sum_i z_i (Wz)_i / sum z^2 with S0 = n_eff (row-standardised)
            Wz = (st["W"] @ Zp.T).T
            return (m / st["nb"]) * (Zp * Wz).sum(1) / zz
        eq = (Zp[:, st["I"]] == Zp[:, st["J"]]) * st["w"]
        return eq.sum(1) / st["nb"]
    dt_ = float if kind == "ord" else int
    obs = {k: stat(z[None, :].astype(dt_), st)[0] for k, st in sets.items()}
    for nullname, bs in [("global", None)] + [(f"block{int(b)}", b) for b in C.BLOCK_SIZES_KM]:
        blocks = None if bs is None else groups_from_ids(grid_ids(d.x_km.values[idx], d.y_km.values[idx], bs))
        use = {k: st for k, st in sets.items() if bs is None or EDGES[k[0] + 1] <= bs}
        nulls = {k: [] for k in use}
        for c0 in range(0, NPERM, CH):
            perm = make_perm(m, blocks, min(CH, NPERM - c0))
            Zp = z[perm].astype(dt_)
            for k, st in use.items():
                nulls[k].append(stat(Zp, st))
        for k, st in use.items():
            nl = np.concatenate(nulls[k]); o = obs[k]; mu = nl.mean(); sd = nl.std(ddof=1)
            row = dict(variable=name, kind=kind, null=nullname, band=f"{EDGES[k[0]]}-{EDGES[k[0]+1]} km", pairs=k[1],
                       n_pairs=st["npairs"], n_samples_with_neighbour=st["nb"], statistic=o, null_mean=mu,
                       null_lo=np.quantile(nl, .025), null_hi=np.quantile(nl, .975), z=(o - mu) / sd,
                       p_perm=(1 + (np.abs(nl - mu) >= abs(o - mu)).sum()) / (len(nl) + 1))
            if kind == "cat":
                row["kappa_s"] = (o - mu) / (1 - mu)
            rows.append(row)
    print(name, round(time.time() - t0), "s", flush=True)
res = pd.DataFrame(rows)
res["q_BH"] = np.nan
for nl_ in res["null"].unique():
    mk = res["null"] == nl_
    res.loc[mk, "q_BH"] = multipletests(res.loc[mk, "p_perm"], method="fdr_bh")[1]
io.save_table(res, f"T03_spatial_autocorrelation{suffix}.csv")
