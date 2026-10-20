"""Run the complete pipeline: python scripts/run_all.py [--sensitivity] [--fast]

Order: 00 data preparation -> 01 contribution -> 02 quality -> 03 spatial autocorrelation ->
04 neighbour consistency -> 05 variance partition -> 06 figures -> session_info.json.
--sensitivity also repeats 03 and 04 with the imputed organic-matter values.
--fast       sets ECHO_FAST=1 (few permutations; for smoke tests only).
"""
import _bootstrap  # noqa: F401
import json, os, platform, subprocess, sys, time
from pathlib import Path

HERE = Path(__file__).resolve().parent
args = sys.argv[1:]
env = dict(os.environ)
if "--fast" in args:
    env["ECHO_FAST"] = "1"
steps = ["00_prepare_data.py", "01_contribution.py", "02_quality.py", "03_spatial_autocorrelation.py",
         "04_neighbour_consistency.py", "05_variance_partition.py", "06_figures.py"]
jobs = [[s] for s in steps]
if "--sensitivity" in args:
    jobs += [["03_spatial_autocorrelation.py", "--om", "imputed"], ["04_neighbour_consistency.py", "--om", "imputed"]]
for job in jobs:
    t0 = time.time()
    print(">>", " ".join(job), flush=True)
    subprocess.run([sys.executable, str(HERE / job[0])] + job[1:], check=True, env=env, cwd=HERE.parent)
    print(f"   done in {time.time()-t0:.0f} s", flush=True)

from echo_ro import config as C, io
import numpy, pandas, scipy, statsmodels, matplotlib
info = dict(python=platform.python_version(), platform=platform.platform(),
            numpy=numpy.__version__, pandas=pandas.__version__, scipy=scipy.__version__,
            statsmodels=statsmodels.__version__, matplotlib=matplotlib.__version__,
            input_sha256=io.sha256(C.INPUT_CSV), seeds=C.SEEDS, fast_mode=C.FAST,
            n_perm=C.N_PERM, n_perm_knn=C.N_PERM_KNN, n_boot=C.N_BOOT, n_boot_knn=C.N_BOOT_KNN)
C.LOGS.mkdir(parents=True, exist_ok=True)
json.dump(info, open(C.LOGS / "session_info.json", "w"), indent=1)
print(json.dumps(info, indent=1))
