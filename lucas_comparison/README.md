# Comparison of ECHO field pH with LUCAS 2018 topsoil pH

`lucas_compare.py` compares cell medians of field pH (ECHO) with laboratory pH in water (LUCAS 2018) on a grid in Stereo 70 (EPSG:3844) anchored at the upper-left corner of the national extent (xmin 134140.3191, ymax 753166.4394), the same grid as in the coverage analysis (`verify_grid.py` checks this).

Design (fixed beforehand): cells with at least 2 ECHO samples and 1 LUCAS point; three pH classes (< 6.0, 6.0–7.5, > 7.5); Spearman ρ of cell medians, linearly weighted kappa, share of cells in the same / an adjacent class, median ECHO − LUCAS difference (descriptive); ρ after centring within blocks of 100 km (separates the regional gradient from local agreement); 95 % intervals by bootstrap over blocks (4,000 replicates). Sensitivity analyses: cells with ≥ 2 collectors; ≥ 3 ECHO samples and ≥ 2 LUCAS points.

```bash
bash lucas_comparison/run_stereo70.sh   # from the repository root; needs data/external/LUCAS-SOIL-2018 RO.csv
```
runs cell sizes 20/30/50 km × blocks 100/200 km (`sens/`), then `origin_robustness.py` (25 shifted grid positions per cell size, `sens_origin/`) and `verify_grid.py`.
A single run: `python lucas_compare.py --echo ../data/input/samples_RO_FINAL.csv --lucas "../data/external/LUCAS-SOIL-2018 RO.csv" --out OUT --cell 30000 --block 100000 --crs 3844 --origin 134140.3191,753166.4394` (see `--help` for options such as `--coords`, `--ph`, `--lat`, `--lon`).

Provided outputs (`sens/*/`): `lucas_summary.json`, `lucas_table.csv`, `lucas_text.txt`, `lucas_fig.*`; combined table `sens/lucas_cell_size_sensitivity.csv`; grid-origin summary `sens_origin/origin_robustness_summary.csv`. The cell-level tables are not distributed because they contain values derived from the restricted LUCAS data.
