#!/usr/bin/env python3
"""Checks that the 30-km Stereo 70 grid used for coverage (RO_grid_count_area.xlsx) and for the LUCAS comparison is the same grid:
cell = floor((x - xmin)/30000), floor((ymax - y)/30000); id = 18*col + row + 1 (extent of the national border in EPSG:3844)."""
import sys, numpy as np, pandas as pd
from pyproj import Transformer
XMIN, YMAX = 134140.3191, 753166.4394
import os
base = sys.argv[1] if len(sys.argv) > 1 else os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "data", "input") + os.sep
g = pd.read_excel(base + "RO_grid_count_area.xlsx"); s = pd.read_excel(base + "RO_voronoi_CLC_altitude.xlsx")
x, y = Transformer.from_crs(4326, 3844, always_xy=True).transform(s.lon.values, s.lat.values)
ids = 18 * np.floor((x - XMIN) / 30000).astype(int) + np.floor((YMAX - y) / 30000).astype(int) + 1
j = pd.concat([pd.Series(ids).value_counts().rename("mine"), g.set_index("id").No_of_points.rename("user")], axis=1).fillna(0)
occ_user, occ_mine = set(g.id[g.No_of_points > 0]), set(pd.Series(ids).unique())
print(f"occupied cells: user {len(occ_user)}, reconstructed {len(occ_mine)}, identical sets: {occ_user == occ_mine}")
print(f"occupied cells with identical counts: {((j.mine == j.user) & (j.user > 0)).sum()} of {(j.user > 0).sum()}; samples in a different cell: {int(abs(j.mine - j.user).sum() / 2)} of {len(s)}")
