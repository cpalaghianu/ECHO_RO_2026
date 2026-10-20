"""Data input/output helpers (no data are modified: derived variables are added as new columns)."""
import json
import hashlib
import numpy as np
import pandas as pd
from . import config as C


def sha256(path):
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for chunk in iter(lambda: f.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def load(om="complete_case"):
    """Load the ECHOREPO export for Romania.

    om = 'complete_case' : soil_organic_matter as exported (6 missing values stay NaN)
    om = 'imputed'       : the 6 missing values replaced by the nearest-3-neighbour value
                           produced by scripts/00_prepare_data.py (sensitivity analysis)
    Adds: dt (UTC timestamp), date (UTC calendar day), month, x_km / y_km (local projection).
    The time-of-day is the upload time and is not used; only the calendar day is analysed.
    """
    d = pd.read_csv(C.INPUT_CSV)
    if om == "imputed":
        imp = pd.read_csv(C.IMPUTED_CSV).set_index("sample_id")["soil_organic_matter"]
        d["soil_organic_matter"] = d["sample_id"].map(imp)
    elif om != "complete_case":
        raise ValueError(om)
    d["dt"] = pd.to_datetime(d["sampling_datetime_utc"], utc=True)
    d["date"] = d["dt"].dt.date
    d["month"] = d["dt"].dt.strftime("%Y-%m")
    lat0 = np.radians(C.LAT0_DEG)
    d["x_km"] = C.EARTH_RADIUS_KM * np.radians(d["lon"].values) * np.cos(lat0)
    d["y_km"] = C.EARTH_RADIUS_KM * np.radians(d["lat"].values)
    return d


def save_json(obj, name):
    C.TABLES.mkdir(parents=True, exist_ok=True)
    with open(C.TABLES / name, "w", encoding="utf-8") as f:
        json.dump(obj, f, indent=1, default=float, ensure_ascii=False)


def save_table(df, name, index=False):
    C.TABLES.mkdir(parents=True, exist_ok=True)
    df.to_csv(C.TABLES / name, index=index)


def rng(module):
    return np.random.default_rng(C.SEEDS[module])
