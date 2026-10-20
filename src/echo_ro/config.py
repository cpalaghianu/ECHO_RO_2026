"""Central configuration: paths, random seeds and analysis parameters.

Every random procedure in the pipeline draws from a generator seeded with the
value stored here, so that all tables and figures can be regenerated exactly
(given the same software versions; see outputs/logs/session_info.json).
Set the environment variable ECHO_FAST=1 for a quick smoke test with reduced
numbers of permutations / bootstrap replicates (results are then approximate).
"""
import os
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
INPUT_CSV = ROOT / "data" / "input" / "samples_RO_FINAL.csv"
DERIVED = ROOT / "data" / "derived"
IMPUTED_CSV = DERIVED / "samples_RO_OMimputed.csv"
IMPUTATION_LOG = DERIVED / "om_imputation_log.csv"
TABLES = ROOT / "outputs" / "tables"
FIGURES = ROOT / "outputs" / "figures"
LOGS = ROOT / "outputs" / "logs"

FAST = os.environ.get("ECHO_FAST", "0") == "1"

# one independent seed per module
SEEDS = {"m01": 20261017, "m02": 20261018, "m03": 20261020,
         "m04": 20261021, "m05": 20261022}

# Organic-matter classes of the ECHO colour chart (% estimate); observed classes only
OM_CLASSES = [1.5, 2.0, 2.5, 3.5, 5.0]

# Coordinates are displaced by a fixed random offset within 1 km (data provider).
# Reported inter-sample distances may therefore be wrong by up to 2 km.
COORD_UNCERTAINTY_KM = 1.0
MIN_PAIR_DISTANCE_KM = 2.0 * COORD_UNCERTAINTY_KM

DIST_EDGES_KM = [2, 5, 10, 20, 50, 100, 200]   # distance bands for spatial statistics
BLOCK_SIZES_KM = [50.0, 100.0]                  # blocks for restricted permutation nulls
CELL_SIZES_KM = [10.0, 25.0]                    # grid cells for ICC / variance partition
KNN_K = [3, 5, 10]
KNN_D0 = [2.0, 5.0]                             # exclusion radii for 'corrected' neighbours

N_PERM = 299 if FAST else 2999                  # permutations (Moran / kappa / operator test)
N_PERM_KNN = 199 if FAST else 1500
N_BOOT = 300 if FAST else 4000                  # bootstrap replicates (proportions, Gini)
N_BOOT_KNN = 200 if FAST else 1500

# Reference latitude for the local equirectangular projection used only to build grid blocks
LAT0_DEG = 46.0
EARTH_RADIUS_KM = 6371.0088
