# Supplementary Material

`make_supplementary.py OUT.docx` builds the Supplementary Material (Tables S1–S7, Figures S1–S3 and the maps) as a Word file from the results in `lucas_comparison/` and `representation_landcover/results/` and from the QGIS maps in `maps_raw/` (JPG). It also writes `lucas_comparison/sens/lucas_cell_size_sensitivity.csv`.

```bash
python supplementary/make_supplementary.py ECHO_RO_Supplementary_Material.docx
```
Requires `python-docx` and `pillow` (`requirements-extra.txt`).
