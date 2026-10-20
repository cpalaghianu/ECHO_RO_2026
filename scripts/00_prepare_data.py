"""Step 0 - integrity check of the input file and derivation of the imputed organic-matter version.

Input : data/input/samples_RO_FINAL.csv  (never modified)
Output: data/derived/samples_RO_OMimputed.csv, data/derived/om_imputation_log.csv

Imputation (used only for the sensitivity analysis and for descriptive completeness):
each missing soil_organic_matter value receives the mean of the three nearest samples with an
observed value (great-circle distance, reported coordinates), rounded to the closest existing
class of the ECHO colour chart (ties resolved towards the lower class). The primary inferential
analyses use complete cases, because neighbour-based imputation builds spatial dependence.
"""
import _bootstrap  # noqa: F401
import numpy as np, pandas as pd
from echo_ro import config as C, io
from echo_ro.stats_utils import haversine_matrix

C.DERIVED.mkdir(parents=True, exist_ok=True)
d = pd.read_csv(C.INPUT_CSV)
print("SHA-256 of input:", io.sha256(C.INPUT_CSV))
assert d.shape == (1091, 11), d.shape
assert not d.sample_id.duplicated().any()

D = haversine_matrix(d.lat.values, d.lon.values)
om0 = d.soil_organic_matter.values.astype(float)
cls = np.array(C.OM_CLASSES)
om = om0.copy()
rows = []
for i in np.where(np.isnan(om0))[0]:
    dist = D[i].copy()
    dist[i] = np.inf
    dist[np.isnan(om0)] = np.inf          # only observed values are used as donors
    nb = np.argsort(dist)[:3]
    m = om0[nb].mean()
    snap = cls[np.argmin(np.abs(cls - m))]  # argmin returns the lower class on ties
    om[i] = snap
    rows.append(dict(sample_id=d.sample_id[i], neighbours=";".join(d.sample_id[nb]),
                     distances_km=";".join(f"{dist[j]:.1f}" for j in nb),
                     neighbour_values=";".join(str(om0[j]) for j in nb),
                     neighbour_mean=round(m, 2), imputed_class=snap))
pd.DataFrame(rows).to_csv(C.IMPUTATION_LOG, index=False)
d2 = d.copy()
d2["soil_organic_matter"] = om
d2["om_imputed"] = np.isnan(om0).astype(int)
d2.to_csv(C.IMPUTED_CSV, index=False)
print(pd.DataFrame(rows).to_string())
