# -*- coding: utf-8 -*-
"""Supplementary Material (Word): Tables S1-S6 and Figures S1-S3 (v1.6: new Table S7 with the reference areas of the altitude bands; Figure S3b recalculated; v1.5: authors, more decimals for values that round to zero, ICC(1) note; v1.4: wording corrections; v1.3: former Figure S4 is now Figure 4 of the article; v1.1: the maps of organic matter, earthworms and pollutants, formerly S5-S7, are now panels b-d of Figure 3 of the article, with the pH map as panel a) for the ECHO Romania article.
Numbering follows the order of first citation in the manuscript. Run: python3 make_supplementary.py OUT.docx
Inputs: lucas_comparison/sens, lucas_comparison/sens_origin, representation_landcover/results, supplementary/maps_raw (QGIS maps, JPG)."""
import sys, json, os, pandas as pd
from PIL import Image, ImageChops
from docx import Document
from docx.shared import Pt, Cm, RGBColor
from docx.enum.text import WD_ALIGN_PARAGRAPH
from docx.enum.section import WD_ORIENT, WD_SECTION
from docx.oxml.ns import qn
from docx.oxml import OxmlElement

OUT = sys.argv[1] if len(sys.argv) > 1 else "ECHO_RO_Supplementary_Material.docx"
HERE = os.path.dirname(os.path.abspath(__file__)); BASE = os.path.join(HERE, "..") + os.sep      # repository root
L = BASE + "lucas_comparison/"; R = BASE + "representation_landcover/results/"; MAPS = BASE + "supplementary/maps_raw/"; CROP = BASE + "supplementary/maps_cropped/"
os.makedirs(CROP, exist_ok=True)
TITLE = "Mapping soil health indicators through citizen science in Romania: a case study of the ECHO network"

def f(x, d=2):
    if x != 0 and d == 2 and round(x, 2) == 0:
        d = 3 if round(x, 3) != 0 else 4      # avoid "-0.000" for values that are small but not zero
    return f"{x:.{d}f}".replace("-", "−")
ci = lambda m: f"{f(m['est'])} ({f(m['lo'])} to {f(m['hi'])})"
def pf(p): return "<0.001" if p < 0.001 else f"{p:.3f}"
def pfl(p): return "p < 0.001" if p < 0.001 else f"p = {p:.3f}"

doc = Document()
st = doc.styles["Normal"]; st.font.name = "Times New Roman"; st.font.size = Pt(10)
st.element.rPr.rFonts.set(qn("w:eastAsia"), "Times New Roman")

def set_section(sec, landscape):
    if landscape:
        sec.orientation = WD_ORIENT.LANDSCAPE; sec.page_width, sec.page_height = Cm(29.7), Cm(21.0)
        sec.left_margin = sec.right_margin = Cm(1.8); sec.top_margin = sec.bottom_margin = Cm(1.6)
    else:
        sec.orientation = WD_ORIENT.PORTRAIT; sec.page_width, sec.page_height = Cm(21.0), Cm(29.7)
        sec.left_margin = sec.right_margin = Cm(2.2); sec.top_margin = sec.bottom_margin = Cm(2.0)

def para(t="", bold=False, size=10, sb=0, sa=4, align=None, italic=False, keep=False):
    p = doc.add_paragraph(); p.paragraph_format.space_before = Pt(sb); p.paragraph_format.space_after = Pt(sa)
    if align: p.alignment = align
    if keep: p.paragraph_format.keep_with_next = True
    if t:
        r = p.add_run(t); r.bold = bold; r.italic = italic; r.font.size = Pt(size)
    return p

def caption(label, title, note=None, keep=False):
    p = para(sb=6, sa=2, keep=keep); r = p.add_run(label + " "); r.bold = True; r.font.size = Pt(10)
    r = p.add_run(title); r.font.size = Pt(10)
    if note:
        n = para(sa=6); rn = n.add_run(note); rn.font.size = Pt(8.5)

def shade(cell, hexcol="E8EEF7"):
    tcPr = cell._tc.get_or_add_tcPr(); sh = OxmlElement("w:shd")
    sh.set(qn("w:val"), "clear"); sh.set(qn("w:color"), "auto"); sh.set(qn("w:fill"), hexcol); tcPr.append(sh)

def table(hdr, data, widths, left_cols=2, fs=8.5):
    t = doc.add_table(rows=1, cols=len(hdr)); t.style = "Table Grid"
    for i, h in enumerate(hdr):
        c = t.rows[0].cells[i]; c.text = ""; r = c.paragraphs[0].add_run(h); r.bold = True; r.font.size = Pt(fs); shade(c)
    trPr = t.rows[0]._tr.get_or_add_trPr(); th = OxmlElement("w:tblHeader"); th.set(qn("w:val"), "true"); trPr.append(th)
    for d in data:
        row = t.add_row().cells
        for i, v in enumerate(d):
            row[i].text = ""; p = row[i].paragraphs[0]; p.paragraph_format.space_after = Pt(0)
            p.alignment = WD_ALIGN_PARAGRAPH.LEFT if i < left_cols else WD_ALIGN_PARAGRAPH.CENTER
            r = p.add_run(str(v)); r.font.size = Pt(fs)
    for row in t.rows:
        for i, w in enumerate(widths): row.cells[i].width = Cm(w)
    return t

def crop_white(src, dst, pad=30, maxw=2600):
    im = Image.open(src).convert("RGB")
    bg = Image.new("RGB", im.size, (255, 255, 255))
    bb = ImageChops.difference(im, bg).point(lambda v: 255 if v > 12 else 0).getbbox()
    x0, y0, x1, y1 = bb; im = im.crop((max(0, x0 - pad), max(0, y0 - pad), min(im.width, x1 + pad), min(im.height, y1 + pad)))
    if im.width > maxw: im = im.resize((maxw, int(im.height * maxw / im.width)), Image.LANCZOS)
    im.save(dst, quality=90); return dst

def picture(path, width_cm, max_h_cm=None):
    im = Image.open(path); h = width_cm * im.height / im.width
    if max_h_cm and h > max_h_cm: width_cm = width_cm * max_h_cm / h
    p = doc.add_paragraph(); p.alignment = WD_ALIGN_PARAGRAPH.CENTER; p.paragraph_format.keep_with_next = True
    p.add_run().add_picture(path, width=Cm(width_cm))

# ------------------------------------------------------------------ cover and list
s0 = doc.sections[0]; set_section(s0, False)
para("Supplementary Material", True, 16, sa=6)
para(TITLE, False, 12, sa=4, italic=True)
para("Palaghianu C., Barnoaiea I., Cosofret C., Savin A.", False, 10, sa=14)
para("Contents", True, 11, sb=6, sa=4)
TABS = [
 ("Table S1", "Land cover and altitude as predictors of the field indicators, with the collector as a random intercept"),
 ("Table S2", "Sensitivity of the agreement between ECHO field pH and LUCAS 2018 pH(H₂O) to cell size, retained cells, block size and position of the grid"),
 ("Table S3", "Median and interquartile range of the field indicators by CORINE Land Cover class and altitude band"),
 ("Table S4", "Pairwise contrasts of the field indicators (Cliff’s δ with collector-cluster bootstrap intervals)"),
 ("Table S5", "Number of ECHO samples in each CORINE Land Cover level-3 class"),
 ("Table S6", "Collector-attributable variance of the field indicators on different scales"),
 ("Table S7", "Representation of the ECHO samples by altitude band, relative to the area of each band in Romania"),
]
FIGS = [
 ("Figure S1", "Sampling locations in the four NUTS 1 macroregions of Romania"),
 ("Figure S2", "Voronoi polygons of the sample locations"),
 ("Figure S3", "Representation of the ECHO samples by CORINE class, altitude band and 30 × 30 km cell"),
]
for a, b in TABS + FIGS:
    p = para(sa=2); r = p.add_run(a + ". "); r.bold = True; p.add_run(b)
para("Sample locations are those stored in ECHOREPO and are displaced by up to 1 km from the true positions; distances are therefore uncertain by about ±2 km. "
     "Maps use the Romanian national boundary and CORINE Land Cover 2018 (vector, version 2020_20u1; European Environment Agency, 2020).", size=9, sb=10, italic=True)

# ------------------------------------------------------------------ tables (landscape)
s1 = doc.add_section(WD_SECTION.NEW_PAGE); set_section(s1, True)

# Table S1: adjusted models
AM = pd.read_csv(R + "O3_adjusted_models.csv"); SM = pd.read_csv(R + "O3_sensitivity_models.csv"); S1m = SM[SM.model.str.startswith("S1")]
nm = {"ph": "Soil pH", "om": "Organic matter", "worms": "Earthworm count", "poll": "Visible pollutants (presence)"}
caption("Table S1.", "Land cover and altitude as predictors of the field indicators, with the collector as a random intercept", keep=True)
rows = []
for _, r in AM.iterrows():
    unit = "OR" if r.indicator == "poll" else "SD"
    if r.indicator == "poll": lrt, s1, s1p = "Wald (GEE)", "–", "–"
    else:
        lrt = (f"χ²(1) = {f(r.LRT_alt_chi2_df1,1)}, {pfl(r.LRT_alt_p)}" if "Altitude" in r.term else f"CORINE class: χ²(2) = {f(r.LRT_clc_chi2_df2,1)}, {pfl(r.LRT_clc_p)}")
        q = S1m[(S1m.indicator == r.indicator) & (S1m.term == r.term)].iloc[0]; s1 = f"{f(q.estimate)} ({f(q.ci_lo)} to {f(q.ci_hi)})"; s1p = pf(q.p_wald)
    rows.append([nm[r.indicator], r.term, f"{int(r.n):,} ({int(r.collectors)})", f"{f(r.estimate)} {unit} ({f(r.ci_lo)} to {f(r.ci_hi)})", pf(r.p_BH), lrt, s1, s1p])
table(["Indicator", "Term (reference: agricultural areas)", "n (collectors)", "Estimate (95% CI)", "p (BH)", "Likelihood-ratio test of the term", "Adding the sampling period: estimate (95% CI)", "p (Wald)"],
      rows, [3.4, 5.0, 2.4, 4.2, 1.6, 4.4, 4.2, 1.6])
para("Linear mixed models (REML) of rank-based normal (Blom) scores, in units of standard deviations (SD), with a random intercept for the collector; the visible-pollutant model is a logistic generalised estimating equation with collector clusters (odds ratio, OR). "
     "Fixed effects: CORINE Land Cover class (artificial surfaces, agricultural areas as reference, forests and semi-natural areas; classes with fewer than 30 samples excluded) and altitude per 100 m. "
     "p-values are Benjamini–Hochberg adjusted over the 12 terms; likelihood-ratio tests used maximum-likelihood fits. The last two columns repeat the models with the sampling period (November–December 2025, January–March, April, May–October 2026) as an additional fixed effect, "
     "which changes the land-cover contrasts but not the altitude effects; the period itself was significant for pH, organic matter and earthworm counts (likelihood-ratio p < 0.00001). "
     "A quadratic altitude term was not significant for pH (p = 0.40) or organic matter (p = 0.64) but was for earthworm counts (p = 0.016).", size=8.5, sb=4)

# Table S2: LUCAS sensitivity
doc.add_page_break()
sets = [("main", "Main (≥2 ECHO, ≥1 LUCAS)"), ("S1_ge2_collectors", "≥2 collectors"), ("S2_ge3echo_ge2lucas", "≥3 ECHO, ≥2 LUCAS")]
rows, csv = [], []
for c in (20, 30, 50):
    a = json.load(open(f"{L}sens/c{c}_b100/lucas_summary.json")); b = json.load(open(f"{L}sens/c{c}_b200/lucas_summary.json"))
    for key, lab in sets:
        m, w = a[key], b[key]
        rows.append([str(c), lab, str(m["n_cells"]), str(m["n_echo_samples"]), str(m["n_lucas_points"]), ci(m["rho"]), ci(m["kappa"]),
                     f"{m['same']['est']:.0f}", f(m["meddiff"]["est"]), ci(m["rho_within"]), ci(w["rho_within"])])
        csv.append(dict(cell_km=c, set=lab, n_cells=m["n_cells"], n_echo=m["n_echo_samples"], n_lucas=m["n_lucas_points"],
                        rho=m["rho"]["est"], rho_lo=m["rho"]["lo"], rho_hi=m["rho"]["hi"], kappa=m["kappa"]["est"], kappa_lo=m["kappa"]["lo"], kappa_hi=m["kappa"]["hi"],
                        same_pct=m["same"]["est"], meddiff=m["meddiff"]["est"], rho_within100=m["rho_within"]["est"], rw100_lo=m["rho_within"]["lo"], rw100_hi=m["rho_within"]["hi"],
                        blocks100=m["n_blocks"], rho_within200=w["rho_within"]["est"], rw200_lo=w["rho_within"]["lo"], rw200_hi=w["rho_within"]["hi"], blocks200=w["n_blocks"]))
pd.DataFrame(csv).to_csv(L + "sens/lucas_cell_size_sensitivity.csv", index=False)
O = pd.read_csv(L + "sens_origin/origin_robustness_summary.csv")
caption("Table S2.", "Sensitivity of the agreement between ECHO field pH and LUCAS 2018 pH(H₂O) to cell size, to the cells retained, to block size and to the position of the grid", keep=True)
para("(a) Cell size, cells retained and block size (Stereo 70 grid, EPSG:3844, anchored at the upper-left corner of the extent of the national boundary: x = 134,140 m, y = 753,166 m)", True, 9, sa=2, keep=True)
table(["Cell (km)", "Cells retained", "n cells", "ECHO", "LUCAS", "ρ (95% CI)", "κw (95% CI)", "Same class (%)", "Median diff.", "ρ within 100-km blocks", "ρ within 200-km blocks"], rows,
      [1.5, 4.0, 1.5, 1.5, 1.6, 3.0, 3.0, 1.9, 1.8, 3.2, 3.2])
para("(b) Position of the grid. The grid was shifted in five steps in each direction over the width of one cell (25 positions per cell size; main analysis; bootstrap reduced to 100 replicates, only point estimates are summarised)", True, 9, sb=8, sa=2, keep=True)
orows = [[str(int(r.cell_km)), f"{int(r.n_cells_min)}–{int(r.n_cells_max)}", f"{f(r.rho_min)} – {f(r.rho_med)} – {f(r.rho_max)}", f"{f(r.kappa_min)} – {f(r.kappa_med)} – {f(r.kappa_max)}", f"{f(r.within_min)} – {f(r.within_med)} – {f(r.within_max)}"] for r in O.itertuples()]
table(["Cell (km)", "n cells (range)", "ρ: minimum – median – maximum", "κw: minimum – median – maximum", "ρ within 100-km blocks: minimum – median – maximum"], orows, [2.0, 3.0, 6.0, 6.0, 8.0])
para("Spearman ρ and linearly weighted kappa (κw, three pH classes: < 6.0, 6.0–7.5, > 7.5) are given with 95% percentile intervals from a bootstrap resampling blocks of 100 km (200 km for the last column of panel a; exploratory, added after the first results). "
     "Within-block ρ is the correlation after centring both cell medians within blocks. Same class: share of cells in the same pH class. Median difference: ECHO − LUCAS (descriptive, not an error estimate). "
     "Blocks of 100 km contain at most four 50-km cells and the 200-km bootstrap rests on 7–10 blocks, so these intervals are approximate. The three cell sizes and the shifted grids use the same samples and are not independent. "
     "The grid of panel (a) is the one used for the coverage analysis (30 km) and was reconstructed from the cell identifiers.", size=8.5, sb=4)

# Table S3: descriptives
doc.add_page_break()
D = pd.read_csv(R + "O3_descriptive_by_group.csv")
caption("Table S3.", "Median and interquartile range of the field indicators by CORINE Land Cover class and altitude band", keep=True)
rows = [[r["factor"], r.level, int(r.n), int(r.collectors), f"{r['pH median']:g} ({r['pH IQR'].replace('-', '–')})", f"{r['OM median']:g} ({r['OM IQR'].replace('-', '–')})",
         f"{r['Earthworms median']:g}", f"{r['Earthworms >=1 (%)']:.0f}", f"{r['Pollutants >=1 (%)']:.0f}", f"{r['altitude median (m)']:.0f}"] for _, r in D.iterrows()]
table(["Factor", "Class or band", "n", "Collectors", "pH: median (IQR)", "Organic matter class: median (IQR)", "Earthworms: median", "Samples with ≥1 earthworm (%)", "Samples with ≥1 visible pollutant (%)", "Median altitude (m)"],
      rows, [3.2, 4.6, 1.4, 1.9, 2.6, 3.2, 2.2, 2.6, 2.8, 2.0])
para("Unweighted descriptive statistics of the 1,091 samples; the samples are clustered within collectors and the figures describe the sample, not the national soil. Organic matter class values are field estimates (1.5, 2, 2.5, 3.5 and 5%). "
     "Wetlands (n = 2) and water bodies (n = 9) are shown for completeness and were not tested. Organic matter was not recorded for six samples, so the number of samples with a value is slightly lower in some rows.", size=8.5, sb=4)

# Table S4: contrasts
P = pd.read_csv(R + "O3_posthoc_contrasts.csv")
ind = {"ph": "pH", "om": "Organic matter class", "worms": "Earthworm count", "poll": "Visible pollutants"}
doc.add_page_break()
caption("Table S4.", "Pairwise contrasts of the field indicators (Cliff’s δ with collector-cluster bootstrap intervals)", keep=True)
rows = [[r.contrast, ind[r.indicator], f"{int(r.n1)} / {int(r.n2)}", f"{r.median1:g} / {r.median2:g}", f"{f(r.cliffs_delta)} ({f(r.ci_lo)} to {f(r.ci_hi)})", pf(r.p_BH_naive),
         "Pre-specified" if r.status == "PRE" else "Exploratory"] for r in P.itertuples()]
table(["Contrast (first vs second group)", "Indicator", "n (first / second)", "Median (first / second)", "Cliff’s δ (95% CI)", "p (BH, naive)", "Status"], rows, [7.0, 3.8, 3.0, 3.0, 4.4, 2.4, 2.6])
para("Cliff’s δ is the probability that a value from the first group exceeds one from the second minus the reverse (negative: lower in the first group); 95% intervals come from a bootstrap that resamples collectors (B = 4,000). "
     "The p-values (Mann–Whitney test, Benjamini–Hochberg adjusted over the 20 contrasts) ignore the clustering of samples within collectors and are therefore anti-conservative; the interpretation rests on the intervals. "
     "Only the contrast of forests with agricultural areas was pre-specified; the others are exploratory.", size=8.5, sb=4)

# Table S5: CLC level 3
C3 = pd.read_csv(R + "O2_representation_clc3_supplementary.csv")
NAMES = {111: "Continuous urban fabric", 112: "Discontinuous urban fabric", 121: "Industrial or commercial units", 122: "Road and rail networks and associated land",
 131: "Mineral extraction sites", 133: "Construction sites", 141: "Green urban areas", 142: "Sport and leisure facilities", 211: "Non-irrigated arable land", 213: "Rice fields",
 221: "Vineyards", 222: "Fruit trees and berry plantations", 231: "Pastures", 242: "Complex cultivation patterns", 243: "Land principally occupied by agriculture, with significant areas of natural vegetation",
 311: "Broad-leaved forest", 312: "Coniferous forest", 313: "Mixed forest", 321: "Natural grasslands", 322: "Moors and heathland", 324: "Transitional woodland-shrub",
 331: "Beaches, dunes, sands", 411: "Inland marshes", 511: "Water courses", 512: "Water bodies"}
doc.add_page_break()
caption("Table S5.", "Number of ECHO samples in each CORINE Land Cover level-3 class", keep=True)
rows = [[int(r["class"]), NAMES[int(r["class"])], int(r.n), f"{r.share_samples_pct:.1f}", r.share_samples_ci.replace(" to ", "–")] for _, r in C3.iterrows()]
table(["Code", "Class", "Samples (n)", "Share of the 1,091 samples (%)", "95% CI (collector bootstrap)"], rows, [1.8, 12.0, 2.6, 4.6, 4.4], left_cols=2)
para("The level-1 classes used in the analysis are the first digit of the code (1 artificial surfaces, 2 agricultural areas, 3 forests and semi-natural areas, 4 wetlands, 5 water bodies). "
     "Level-3 classes without samples are not listed. The classes have no reference areas at this level, so only the number of samples is reported.", size=8.5, sb=4)

# Table S6: robustness of the variance shares to the scale of the indicator
T9 = pd.read_csv(R + "T09_variance_share_by_scale.csv")
doc.add_page_break()
caption("Table S6.", "Collector-attributable variance of the field indicators on different scales", keep=True)
pc = lambda x: f"{100*x:.1f}"
rows = [[r.variable, r.scale, int(r.n), pc(r.share_collector), pc(r.share_cell25), pc(r.share_block100), pc(r.share_residual), pf(r.p_collector), f(r.ICC1_collector)] for r in T9.itertuples()]
table(["Indicator", "Scale of the response", "n", "Collector (%)", "25-km cell (%)", "100 × 100 km block (%)", "Residual (%)", "p (collector)", "ICC(1), collector"],
      rows, [3.6, 5.4, 1.4, 2.4, 2.4, 3.0, 2.4, 2.4, 2.6], left_cols=2)
para("Crossed random-effects model of the main analysis (collector + 25-km cell + 100 × 100 km block + residual; REML) refitted to the rank-based normal (Blom) scores used in the article (reference), to mid-ranks expressed as z-scores, "
     "to the original values (standardised) and, for earthworm counts, to log(1 + count) and the square root of the count. Shares are ratios of variance components and can be compared across scales. "
     "p: likelihood-ratio test of the collector component with the 50:50 mixture of χ² distributions. ICC(1): one-way intraclass correlation among the 212 collectors with at least two samples (859 samples), computed on the scale shown; the article reports it on mid-ranks (0.47, 0.21 and 0.34), which is why the value for the Blom scores is slightly higher. The collector share ranged from 40 to 47% for pH, from 13 to 18% for organic matter and from 30 to 31% for earthworm counts; "
     "the ordering of the indicators did not change. Descriptive: no intervals are given, and the Gaussian models remain an approximation for ordinal classes; ordinal mixed models were not fitted.", size=8.5, sb=4)

# Table S7: representation by altitude band (reference areas from the digital elevation model)
RA = pd.read_csv(R + "O2_representation_altitude.csv")
REFA = pd.read_csv(R + "../ref_alt.csv").set_index("band").area_ha
doc.add_page_break()
caption("Table S7.", "Representation of the ECHO samples by altitude band, relative to the area of each band in Romania", keep=True)
lab = {"<200": "< 200", "200-499": "200–499", "500-799": "500–799", "800-1199": "800–1,199", ">=1200": "≥ 1,200"}
rows = [[lab[r["class"]], int(r.n), f"{r.share_samples_pct:.1f}", f"{r.voronoi_share_pct:.1f}", f"{REFA[r['class']]/100:,.0f}", f"{r.ref_share_pct:.1f}",
         f"{r.RR_vs_reference:.2f}", r.RR_ci.replace(" to ", "–")] for _, r in RA.iterrows()]
table(["Altitude band (m)", "n", "Samples (%)", "Voronoi area (%)", "Area of the band in Romania (km²)", "Area of Romania (%)", "Representation ratio", "95% CI"],
      rows, [3.4, 1.6, 2.4, 2.8, 4.4, 3.0, 3.0, 3.0], left_cols=1)
para("Representation ratio = share of the samples / share of the area of Romania in the band. The area of each band was calculated from the Copernicus DEM GLO-30 (total 238,394 km²); "
     "intervals are 95% percentile intervals of a bootstrap that resamples collectors (4,000 replicates). The Voronoi share is the area of the polygons of the samples in the band, used as a check on clustering.", size=8.5, sb=4)

# ------------------------------------------------------------------ figures (portrait)
s2 = doc.add_section(WD_SECTION.NEW_PAGE); set_section(s2, False)
maps = {"S1": ("nuts.jpg", "Sampling locations of the 1,091 ECHO soil samples within the four NUTS 1 macroregions of Romania (RO1–RO4). Locations are those stored in ECHOREPO, displaced by up to 1 km."),
        "S2": ("voronoi.jpg", "Voronoi (Thiessen) polygons of the 1,091 sample locations, clipped to the national boundary. Each polygon contains the points that are closer to its sample than to any other and is used as the footprint of the sample in the analysis of representation (median footprint 134 km²; Gini coefficient of the footprints 0.57; Kish effective sample size 372, i.e. 34% of the samples).")}
def fig_map(code, first=False):
    src, cap = maps[code]; dst = CROP + src.replace(" ", "_")
    crop_white(MAPS + src, dst)
    if not first: doc.add_page_break()
    caption(f"Figure {code}.", cap, keep=True); picture(dst, 16.5, 19.0)

fig_map("S1", True)
fig_map("S2")
# S3 on its own page
doc.add_page_break()
caption("Figure S3.", "Representation of the ECHO samples: (a) share of the samples, of the Voronoi area and of the true area (CORINE Land Cover 2018) by CORINE level-1 class (error bars: 95% intervals of a bootstrap resampling collectors; wetlands n = 2 and water bodies n = 9); "
        "(b) share of the samples, of the Voronoi area and of the true area (from the digital elevation model) by altitude band; (c) number of 30 × 30 km cells by number of samples, observed and expected if the same number of samples were placed at random (mean of 10,000 simulations).", keep=True)
picture(R + "Fig_O2_representation.png", 16.5)
doc.core_properties.title = "Supplementary Material – " + TITLE; doc.core_properties.author = "Palaghianu C. et al."; doc.core_properties.comments = ""; import datetime as _dt; doc.core_properties.created = doc.core_properties.modified = _dt.datetime(2026, 10, 20)
doc.save(OUT); print("saved", OUT)
