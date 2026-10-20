# Input data

| File | Content | In this repository? |
|---|---|---|
| `samples_RO_FINAL.csv` | Romanian records exported from ECHOREPO (1,091 samples × 11 columns; coordinates displaced by the data provider by a fixed random offset ≤ 1 km; collectors appear as pseudonymous identifiers). Used by `scripts/` and `lucas_comparison/`. | see note below |
| `RO_voronoi_CLC_altitude.xlsx` | One row per sample with indicators, collector, Voronoi-polygon area (ha), CORINE Land Cover 2018 code and altitude (m); built by the authors in GIS from the file above. Used by `representation_landcover/` and `lucas_comparison/verify_grid.py`. | see note below |
| `RO_grid_count_area.xlsx` | 30 × 30 km Stereo 70 (EPSG:3844) grid cells intersecting Romania: number of samples and area (ha); last row = totals. | yes |
| `nuts_romania_count.xlsx` | Number of samples per NUTS 2 region and NUTS 3 county, with the Eurostat typology codes. | yes |
| `procente_altitudine.xlsx` | Area (ha) of the five altitude bands (< 200, 200–499, 500–799, 800–1,199, ≥ 1,200 m) of Romania. | yes |

The two sample-level files are *not* part of the repository unless the authors added them (check this folder). They can be rebuilt as follows:

* **ECHOREPO export.** Download the ECHOREPO dataset (Osychenko, 2026, Zenodo, https://doi.org/10.5281/zenodo.19722513), keep the Romanian records, and save them with the 11 columns listed above as `samples_RO_FINAL.csv`. The file analysed in the article has SHA-256 `f421dd917c506b06173cf5f9ddb020603f390e35e01ed41e1cc1132cf9cde2e6` (`scripts/run_all.py` writes the hash of the file in use to `outputs/logs/session_info.json`).
* **Voronoi/CORINE/altitude table.** Voronoi polygons of the sample points clipped to the national boundary, intersected with CORINE Land Cover 2018 (level 3 code in `Code_18`, level 1 class in `Corine2018`) and a digital elevation model; one row per sample, area in `voronoi_ar`.

Reference areas used for the representation ratios: `representation_landcover/ref_clc.csv` (CORINE level-1 classes) and `representation_landcover/ref_alt.csv` (altitude bands, from `procente_altitudine.xlsx`).
