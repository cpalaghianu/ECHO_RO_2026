# echo_RO – analysis code for the Romanian ECHO citizen-science soil data

Code, derived inputs, result tables and figures for the article

> Palaghianu C., Barnoaiea I., Cosofret C., Savin A. *Mapping soil health indicators through citizen science in Romania: a case study of the ECHO network.* AgroLife Scientific Journal (submitted).

The analysis uses the Romanian records of ECHOREPO (1,091 samples from 444 collectors) and addresses contribution and data quality, spatial autocorrelation, collector effects, spatial representation (grid, Voronoi polygons, NUTS units, CORINE land cover, altitude) and the comparison of field pH with the LUCAS 2018 topsoil survey.

**Archived release:** [![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.23496521.svg)](https://doi.org/10.5281/zenodo.23496521) **Repository:** https://github.com/cpalaghianu/ECHO_RO_2026

## Repository layout

| Path | Content |
|---|---|
| `scripts/` | main pipeline, steps `00`–`09` (data checks, contribution, quality, Moran's I, neighbour consistency, variance partition, figures, sensitivity, NUTS representation, robustness to scale); `run_all.py` runs steps 0–6 |
| `src/echo_ro/` | configuration (seeds, parameters), data loader, statistical helpers, figure style |
| `tests/check_reference_values.py` | compares key results with the reference run |
| `outputs/` | result tables (`tables/`), figures (`figures/`, PDF and 300-dpi PNG) and run logs (`logs/`) |
| `representation_landcover/` | spatial representation, indicators by land cover and altitude; own tests (`test_o2_o3.py`) and results (`results/`) |
| `lucas_comparison/` | comparison of field pH with LUCAS 2018 on a Stereo 70 grid; cell-size, block-size and grid-origin sensitivity (`sens/`, `sens_origin/`) |
| `supplementary/` | generator of the Supplementary Material (Word) and the QGIS maps it uses |
| `data/` | input data (see `data/input/README.md`), derived files, external data (see `data/external/README.md`) |
| `docs/` | proposed Methods text and a statistical review of the methods |

## Requirements and quick start

Python 3.13; package versions used for the reported results are pinned in `requirements.txt` (core pipeline) and `requirements-extra.txt` (representation, LUCAS and supplementary scripts).

```bash
python -m venv .venv && source .venv/bin/activate      # optional
pip install -r requirements.txt -r requirements-extra.txt
# put the input files in data/input/ (see data/input/README.md)

python scripts/run_all.py                 # steps 0-6, ~12 min on 2 CPU cores
python scripts/run_all.py --sensitivity   # also repeats steps 3-4 with imputed organic matter
python scripts/run_all.py --fast          # smoke test with few permutations (approximate)
python scripts/07_sensitivity.py          # further sensitivity analyses
python scripts/08_nuts_representation.py  # NUTS 2 / NUTS 3 representation
python scripts/09_robustness_scale.py     # variance shares on alternative scales
python tests/check_reference_values.py    # compare with reference values

bash representation_landcover/run.sh      # representation, land cover, altitude (+ 40 consistency checks)
bash lucas_comparison/run_stereo70.sh     # LUCAS comparison (needs data/external/, ~15 min)
python supplementary/make_supplementary.py ECHO_RO_Supplementary_Material.docx
```

Each script can also be run alone from the repository root; they need only the outputs of earlier steps.

## Reproducibility

* All random procedures use seeded NumPy generators (seeds in `src/echo_ro/config.py`); results are identical on re-runs with the same package versions. Other versions may change Monte-Carlo p-values in the third decimal.
* Parameters (distance bands, numbers of permutations and bootstrap replicates, block sizes) are in `src/echo_ro/config.py`.
* Coordinates in ECHOREPO are displaced by the data provider by a fixed random offset ≤ 1 km; pairs closer than 2 km are excluded from spatial statistics and distances are uncertain by up to 2 km.
* `outputs/logs/session_info.json` records software versions, seeds and the SHA-256 of the input file.

## Data

* ECHOREPO: Osychenko, O. (2026). *ECHOREPO Soil Dataset: Soil observations, parameters, images and biodiversity* [Data set]. Zenodo. https://doi.org/10.5281/zenodo.19722513
* LUCAS 2018 Topsoil: European Soil Data Centre (ESDAC), access on request; **not redistributed** here.
* CORINE Land Cover 2018 and the digital elevation model were processed in GIS by the authors; the derived tables are listed in `data/input/README.md`.

## Citation

Please cite the article and the archived release of this repository (see `CITATION.cff`).

## Licence

See the `LICENSE` file.
