#!/bin/bash
# Spatial representation (grid, Voronoi, CORINE, altitude) and indicators vs land cover / altitude.
# Run from the repository root:  bash representation_landcover/run.sh
# Needs data/input/RO_voronoi_CLC_altitude.xlsx (sample-level table with Voronoi area, CORINE code and altitude; see data/input/README.md).
set -e
cd "$(dirname "$0")"
python3 echo_o2_o3_analysis.py --samples ../data/input/RO_voronoi_CLC_altitude.xlsx --grid ../data/input/RO_grid_count_area.xlsx \
        --ref-clc ref_clc.csv --ref-alt ref_alt.csv --out results
python3 test_o2_o3.py
