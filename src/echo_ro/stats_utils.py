"""Statistical helper functions used by several analysis scripts."""
import numpy as np
import pandas as pd
from scipy import stats
from . import config as C


def haversine_matrix(lat, lon):
    """Great-circle distance matrix (km) between all pairs of points."""
    la, lo = np.radians(lat), np.radians(lon)
    dlat = la[:, None] - la[None, :]
    dlon = lo[:, None] - lo[None, :]
    a = np.sin(dlat / 2) ** 2 + np.cos(la)[:, None] * np.cos(la)[None, :] * np.sin(dlon / 2) ** 2
    return 2 * C.EARTH_RADIUS_KM * np.arcsin(np.sqrt(a))


def grid_ids(x_km, y_km, size_km):
    """Integer id of the square grid cell (size_km) containing each point."""
    ix = np.floor(x_km / size_km).astype(np.int64)
    iy = np.floor(y_km / size_km).astype(np.int64)
    return pd.factorize(ix * 1_000_003 + iy)[0]


def groups_from_ids(ids):
    return [np.where(ids == g)[0] for g in np.unique(ids)]


def gini(x):
    """Gini coefficient of non-negative values (0 = equal, 1 = maximal concentration)."""
    x = np.sort(np.asarray(x, float))
    k = len(x)
    cum = np.cumsum(x)
    return (k + 1 - 2 * cum.sum() / cum[-1]) / k


def hill_numbers(counts):
    """Hill numbers of order 0, 1, 2 (effective number of contributors)."""
    p = np.asarray(counts, float) / np.sum(counts)
    return {"N0": len(p), "N1": float(np.exp(-(p * np.log(p)).sum())), "N2": float(1.0 / (p ** 2).sum())}


def cluster_bootstrap_ratio(num, den, rng, B):
    """Percentile CI (95 %) for sum(num)/sum(den) resampling clusters (collectors) with replacement."""
    num = np.asarray(num, float)
    den = np.asarray(den, float)
    g = len(num)
    out = np.empty(B)
    for b in range(B):
        s = rng.integers(0, g, g)
        out[b] = num[s].sum() / den[s].sum()
    return np.quantile(out, [0.025, 0.975])


def blom_inverse_normal(x):
    """Rank-based inverse-normal transform (Blom offsets); ties receive the mean rank."""
    r = stats.rankdata(x)
    return stats.norm.ppf((r - 3.0 / 8.0) / (len(x) + 1.0 / 4.0))


def icc1(rank, group):
    """One-way random-effects ICC(1) for unbalanced groups (groups with a single value are dropped)."""
    g = pd.Series(group)
    keep = (g.map(g.value_counts()) >= 2).values
    r = np.asarray(rank, float)[keep]
    gg = pd.factorize(g[keep])[0]
    k = gg.max() + 1
    n = len(r)
    if k < 2:
        return np.nan
    sums = np.bincount(gg, r)
    ns = np.bincount(gg).astype(float)
    gm = r.mean()
    ssb = (ns * (sums / ns - gm) ** 2).sum()
    ssw = ((r - (sums / ns)[gg]) ** 2).sum()
    msb, msw = ssb / (k - 1), ssw / (n - k)
    k0 = (n - (ns ** 2).sum() / n) / (k - 1)
    return (msb - msw) / (msb + (k0 - 1) * msw)
