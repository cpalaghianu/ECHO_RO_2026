# Spatial representation, land cover and altitude

`echo_o2_o3_analysis.py` — (i) coverage of Romania by the 30 × 30 km Stereo 70 grid and by Voronoi polygons of the sample points; (ii) representation ratios (share of samples or of Voronoi area ÷ share of national area) by CORINE Land Cover level-1 class and by altitude band, with cluster (collector) bootstrap intervals; (iii) field indicators by land-cover class and altitude band: Kruskal–Wallis ε², Cliff's δ, collector-adjusted models (mixed models, GEE) and sensitivity models. Collectors are always treated as clusters.

```bash
bash representation_landcover/run.sh        # from the repository root
```
runs the analysis with the reference areas `ref_clc.csv` and `ref_alt.csv` and then `test_o2_o3.py` (40 checks: independent re-computation with pandas/scipy of the main quantities, consistency of the reference areas). Results go to `results/` (CSV, Excel workbook, figures).

| File | Content |
|---|---|
| `ref_clc.csv` | CORINE Land Cover 2018 level-1 areas of Romania (ha), computed by the authors in GIS |
| `ref_alt.csv` | areas (ha) of the five altitude bands, from `data/input/procente_altitudine.xlsx` (total 23,839,399 ha) |
| `lucas_landcover_context.py` | the same land-cover contrasts in LUCAS 2018 (descriptive context; needs `data/external/LUCAS-SOIL-2018 RO.csv`): `python lucas_landcover_context.py "../data/external/LUCAS-SOIL-2018 RO.csv" results` |

Inputs: `data/input/RO_voronoi_CLC_altitude.xlsx` and `RO_grid_count_area.xlsx` (see `data/input/README.md`).
