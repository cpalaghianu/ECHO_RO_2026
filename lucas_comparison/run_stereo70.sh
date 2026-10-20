#!/bin/bash
# ECHO-LUCAS comparison on the Stereo 70 grid (EPSG:3844) anchored at the upper-left corner of the national extent
# (xmin 134140.3191, ymax 753166.4394). Six runs (cell 20/30/50 km x block 100/200 km) into sens/, then the 25 grid
# positions per cell size (origin_robustness.py -> sens_origin/), then a check of the coverage grid (verify_grid.py).
# Run from the repository root:  bash lucas_comparison/run_stereo70.sh
# Needs data/external/LUCAS-SOIL-2018 RO.csv (LUCAS 2018 topsoil, Romania; request access from ESDAC - not redistributed here).
set -e
cd "$(dirname "$0")"
E="../data/input/samples_RO_FINAL.csv"; L="../data/external/LUCAS-SOIL-2018 RO.csv"
mkdir -p sens
for c in 20 30 50; do for b in 100 200; do
  python3 lucas_compare.py --echo "$E" --lucas "$L" --out sens/c${c}_b${b} --cell ${c}000 --block ${b}000 --crs 3844 --origin 134140.3191,753166.4394 > sens/log_c${c}_b${b}.txt 2>&1 &
done; done
wait
python3 origin_robustness.py "$E" "$L" > origin_robustness.log 2>&1
python3 verify_grid.py     # needs data/input/RO_voronoi_CLC_altitude.xlsx and RO_grid_count_area.xlsx
