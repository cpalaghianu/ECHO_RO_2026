#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
ECHO (Romania, v11) vs LUCAS Topsoil 2018: pH(H2O) at the level of 30 x 30 km grid cells.

Usage (minimal):
  python3 lucas_compare.py --echo ../data/input/samples_RO_FINAL.csv --lucas LUCAS-SOIL-2018.csv --out out_lucas
Optional:
  --coords LUCAS_coordinates.csv   separate coordinate file (joined on POINTID)
  --lat/--lon/--ph/--country       force column names if auto-detection fails
  --simulated                      stamps figures "SIMULATED DATA" (use for tests only)

What it does (all choices fixed in advance):
  1. LUCAS: keep Romania, valid pH(H2O); ECHO: field pH 'soil_ph_field'.
  2. Both datasets -> EPSG:3844 (Stereo 70) -> 30 km cells anchored at the national extent (xmin 134140.3191, ymax 753166.4394); same grid as the coverage analysis.
  3. Cell medians; keep cells with >=2 ECHO samples and >=1 LUCAS point.
  4. Spearman rho on cell medians; weighted (linear) kappa on 3 pH classes (<6.0; 6.0-7.5; >7.5);
     % same class / within one class; median difference ECHO - LUCAS (descriptive only);
     rho after centring both medians within 100-km blocks (removes the regional gradient).
  5. 95% CIs: percentile bootstrap resampling 100-km blocks.
  6. Sensitivities: (S1) cells with >=2 different ECHO collectors; (S2) >=3 ECHO and >=2 LUCAS points.
"""
import argparse, json, os, re, sys, warnings
import numpy as np, pandas as pd
from scipy.stats import spearmanr
from pyproj import Transformer

ROMANIA_AREA_KM2 = 238397.0
CLASS_BREAKS = (6.0, 7.5)           # <6.0 | 6.0-7.5 | >7.5
CLASS_LABELS = ['acid (<6.0)', 'slightly acid-neutral (6.0-7.5)', 'alkaline (>7.5)']

def die(msg):
    print('ERROR: ' + msg, file=sys.stderr); sys.exit(1)

def read_any(path):
    for enc in ('utf-8-sig', 'latin-1'):
        for sep in (',', ';', '\t'):
            try:
                d = pd.read_csv(path, sep=sep, encoding=enc, low_memory=False)
                if d.shape[1] > 2:
                    return d
            except Exception:
                pass
    die(f'cannot read {path}')

def pick(cols, names):
    low = {c.lower().strip(): c for c in cols}
    for n in names:
        if n in low:
            return low[n]
    return None

def detect(lucas, args):
    cols = list(lucas.columns)
    ph = args.ph or pick(cols, ['ph_h2o', 'ph(h2o)', 'ph h2o', 'ph_in_h2o', 'ph_water'])
    if ph is None:
        c = [x for x in cols if 'ph' in x.lower() and 'h2o' in x.lower().replace('₂', '2')]
        ph = c[0] if c else None
    lat = args.lat; lon = args.lon
    if not (lat and lon):
        for a, b in [('gps_lat', 'gps_long'), ('gps_lat', 'gps_lon'), ('th_lat', 'th_long'), ('th_lat', 'th_lon'),
                     ('lat', 'lon'), ('lat', 'long'), ('latitude', 'longitude')]:
            la, lo = pick(cols, [a]), pick(cols, [b])
            if la and lo:
                lat, lon = la, lo; break
    ctry = args.country or pick(cols, ['nuts0', 'nuts_0', 'cntr_code', 'cntr', 'country', 'survey_country', 'ctry'])
    pid = pick(cols, ['pointid', 'point_id', 'id'])
    return ph, lat, lon, ctry, pid

def to_num(s):
    return pd.to_numeric(s.astype(str).str.replace(',', '.', regex=False).str.strip(), errors='coerce')

GRID = {'crs': 3844, 'origin': (134140.3191, 753166.4394)}      # set in main() from --crs / --origin

def cell_ids(lat, lon, cell):
    t = Transformer.from_crs(4326, GRID['crs'], always_xy=True)
    x, y = t.transform(np.asarray(lon, float), np.asarray(lat, float))
    ix = np.floor((x - GRID['origin'][0]) / cell).astype(int); iy = np.floor((y - GRID['origin'][1]) / cell).astype(int)
    return ix, iy

def classify(v):
    v = np.asarray(v, float)
    return np.where(v < CLASS_BREAKS[0], 0, np.where(v <= CLASS_BREAKS[1], 1, 2))

def kappa_lin(a, b, k=3):
    O = np.zeros((k, k))
    np.add.at(O, (a, b), 1)
    n = O.sum()
    if n == 0: return np.nan
    E = np.outer(O.sum(1), O.sum(0)) / n
    w = np.abs(np.subtract.outer(np.arange(k), np.arange(k))) / (k - 1)
    den = (w * E).sum()
    return np.nan if den == 0 else 1 - (w * O).sum() / den

def safe_rho(a, b):
    if len(a) < 4 or np.ptp(a) == 0 or np.ptp(b) == 0: return np.nan
    return spearmanr(a, b)[0]

def metrics(e, l, ec, lc):
    ca, cb = classify(e), classify(l)
    return dict(rho=safe_rho(e, l), kappa=kappa_lin(ca, cb),
                same=np.mean(ca == cb) * 100, within1=np.mean(np.abs(ca - cb) <= 1) * 100,
                meddiff=np.median(e - l), rho_within=safe_rho(ec, lc))

def analyse(df, boot, rng):
    df = df.reset_index(drop=True)
    e, l = df.ph_e.values, df.ph_l.values
    # centre within 100-km blocks (blocks with one cell contribute 0 and are dropped)
    gm_e = df.groupby('blk').ph_e.transform('mean').values; gm_l = df.groupby('blk').ph_l.transform('mean').values
    gn = df.groupby('blk').ph_e.transform('size').values
    ec, lc = (e - gm_e), (l - gm_l)
    keep = gn >= 2
    obs = metrics(e, l, ec[keep], lc[keep])
    blocks = [np.where(df.blk.values == b)[0] for b in df.blk.unique()]
    nb = len(blocks)
    names = ['rho', 'kappa', 'same', 'within1', 'meddiff', 'rho_within']
    draws = {k: [] for k in names}
    with warnings.catch_warnings():
        warnings.simplefilter('ignore')
        for _ in range(boot):
            idx = np.concatenate([blocks[i] for i in rng.integers(0, nb, nb)])
            k2 = keep[idx]
            m = metrics(e[idx], l[idx], ec[idx][k2], lc[idx][k2])
            for k in names: draws[k].append(m[k])
    out = dict(n_cells=int(len(df)), n_blocks=int(nb), n_echo_samples=int(df.n_e.sum()), n_lucas_points=int(df.n_l.sum()),
               echo_median=float(np.median(e)), lucas_median=float(np.median(l)))
    for k in names:
        d = np.array(draws[k], float); d = d[~np.isnan(d)]
        out[k] = dict(est=None if np.isnan(obs[k]) else float(obs[k]),
                      lo=float(np.percentile(d, 2.5)) if len(d) > 50 else None,
                      hi=float(np.percentile(d, 97.5)) if len(d) > 50 else None)
    cm = np.zeros((3, 3), int); np.add.at(cm, (classify(e), classify(l)), 1)
    out['confusion_echo_rows_lucas_cols'] = cm.tolist()
    return out

def fmt(x, d=2):
    return 'NA' if x is None else f'{x:.{d}f}'

def ci(m, d=2):
    return f"{fmt(m['est'], d)} (95% CI {fmt(m['lo'], d)} to {fmt(m['hi'], d)})"

def make_figure(df, res, outdir, simulated):
    import matplotlib; matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    plt.rcParams.update({'font.family': 'serif', 'font.serif': ['Liberation Serif', 'Times New Roman', 'DejaVu Serif'],
                         'font.size': 9, 'axes.spines.top': False, 'axes.spines.right': False})
    BLUE, INK, MUTED = '#2a78d6', '#0b0b0b', '#52514e'
    fig, ax = plt.subplots(1, 2, figsize=(7.2, 3.3), gridspec_kw={'width_ratios': [1.15, 1]})
    rng = np.random.default_rng(1)
    a = ax[0]
    jit = rng.uniform(-0.12, 0.12, len(df))
    a.scatter(df.ph_l, df.ph_e + jit, s=14 + 4 * np.clip(df.n_e, 0, 10), color=BLUE, alpha=0.75,
              edgecolor='white', linewidth=0.5, zorder=3)
    lim = [min(df.ph_l.min(), df.ph_e.min()) - 0.3, max(df.ph_l.max(), df.ph_e.max()) + 0.3]
    a.plot(lim, lim, color=MUTED, lw=0.8, ls='--', zorder=1)
    a.set_xlim(lim); a.set_ylim(lim)
    a.set_xlabel('LUCAS 2018 pH (H$_2$O), cell median'); a.set_ylabel('ECHO field pH, cell median')
    r = res['main']
    a.text(0.03, 0.97, f"n = {r['n_cells']} cells\n$\\rho$ = {fmt(r['rho']['est'])} ({fmt(r['rho']['lo'])} to {fmt(r['rho']['hi'])})",
           transform=a.transAxes, va='top', ha='left', color=INK, fontsize=8.5)
    a.set_title('(a)', loc='left', fontsize=9, color=INK)
    b = ax[1]
    cm = np.array(r['confusion_echo_rows_lucas_cols'])
    b.imshow(cm, cmap='Blues', vmin=0, vmax=max(cm.max(), 1))
    for i in range(3):
        for j in range(3):
            b.text(j, i, str(cm[i, j]), ha='center', va='center',
                   color='white' if cm[i, j] > 0.55 * max(cm.max(), 1) else INK, fontsize=9)
    b.set_xticks(range(3)); b.set_yticks(range(3))
    b.set_xticklabels(['<6.0', '6.0–7.5', '>7.5']); b.set_yticklabels(['<6.0', '6.0–\n7.5', '>7.5'])
    b.set_xlabel('LUCAS class'); b.set_ylabel('ECHO class', labelpad=8)
    b.set_title(f"(b)  weighted $\\kappa$ = {fmt(r['kappa']['est'])}", loc='left', fontsize=9, color=INK)
    for s in b.spines.values(): s.set_visible(False)
    if simulated:
        fig.text(0.5, 0.5, 'SIMULATED DATA – NOT FOR PUBLICATION', color='red', alpha=0.35, fontsize=22,
                 rotation=20, ha='center', va='center')
    fig.tight_layout()
    for ext in ('png', 'pdf'):
        fig.savefig(os.path.join(outdir, f'lucas_fig.{ext}'), dpi=300)
    plt.close(fig)

def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('--echo', required=True); ap.add_argument('--lucas', required=True)
    ap.add_argument('--coords'); ap.add_argument('--out', default='out_lucas')
    ap.add_argument('--cell', type=float, default=30000.0); ap.add_argument('--boot', type=int, default=4000)
    ap.add_argument('--seed', type=int, default=20261018); ap.add_argument('--block', type=float, default=100000.0, help='block size (m) for centring and bootstrap; default 100 km')
    ap.add_argument('--min-echo', type=int, default=2); ap.add_argument('--min-lucas', type=int, default=1)
    ap.add_argument('--lat'); ap.add_argument('--lon'); ap.add_argument('--ph'); ap.add_argument('--country')
    ap.add_argument('--simulated', action='store_true')
    ap.add_argument('--crs', type=int, default=3844, help='EPSG code of the grid (default 3844 = Stereo 70; 3035 = ETRS89-LAEA)')
    ap.add_argument('--origin', default='134140.3191,753166.4394', help='x,y (m) of a grid node in the grid CRS; default = xmin, ymax of the national extent (Stereo 70)')
    args = ap.parse_args()
    GRID['crs'] = args.crs; GRID['origin'] = tuple(float(v) for v in args.origin.split(','))
    os.makedirs(args.out, exist_ok=True)
    rng = np.random.default_rng(args.seed)

    # ---------------- ECHO
    echo = read_any(args.echo)
    need = {'lat', 'lon', 'soil_ph_field', 'collected_by'}
    if not need <= set(echo.columns): die(f'ECHO file lacks columns {need - set(echo.columns)}')
    echo = echo.dropna(subset=['lat', 'lon', 'soil_ph_field']).copy()
    echo['ix'], echo['iy'] = cell_ids(echo.lat, echo.lon, args.cell)
    echo['cell'] = echo.ix.astype(str) + '_' + echo.iy.astype(str)

    # ---------------- LUCAS
    lucas = read_any(args.lucas)
    ph, lat, lon, ctry, pid = detect(lucas, args)
    if args.coords:
        co = read_any(args.coords)
        cph, clat, clon, _, cpid = detect(co, args)
        if not (clat and clon and cpid and pid): die(f'coords file: cannot find id/lat/lon. Columns: {list(co.columns)}')
        lucas = lucas.merge(co[[cpid, clat, clon]].rename(columns={cpid: pid}), on=pid, how='left', suffixes=('', '_c'))
        lat = clat + '_c' if clat in lucas.columns and clat + '_c' in lucas.columns else clat
        lon = clon + '_c' if clon in lucas.columns and clon + '_c' in lucas.columns else clon
    print(f'LUCAS columns used -> pH: {ph} | lat: {lat} | lon: {lon} | country: {ctry} | id: {pid}')
    if ph is None: die(f'pH(H2O) column not found. Use --ph. Columns: {list(lucas.columns)}')
    if not (lat and lon): die(f'coordinate columns not found. Use --lat/--lon or --coords. Columns: {list(lucas.columns)}')
    n_all = len(lucas)
    depth_counts = lucas[pick(list(lucas.columns), ['depth'])].value_counts().to_dict() if pick(list(lucas.columns), ['depth']) else {}
    if ctry:
        lucas = lucas[lucas[ctry].astype(str).str.strip().str.upper().isin(['RO', 'ROU', 'ROMANIA'])]
    else:
        print('WARNING: no country column; using a bounding box (includes some foreign points!)')
    lucas = lucas.copy()
    lucas['ph'] = to_num(lucas[ph]); lucas['la'] = to_num(lucas[lat]); lucas['lo'] = to_num(lucas[lon])
    n_ro = len(lucas)
    bad_ph = int(lucas.ph.isna().sum())
    lucas = lucas.dropna(subset=['ph', 'la', 'lo'])
    if not ((lucas.la.between(-90, 90)).all() and (lucas.lo.between(-180, 180)).all()):
        die('coordinates are not in degrees (WGS84). Convert first or tell me the CRS.')
    if not ctry:
        lucas = lucas[lucas.la.between(43.5, 48.4) & lucas.lo.between(20.2, 29.8)]
    if not lucas.ph.between(2.5, 10.5).all():
        print('WARNING: some pH values outside 2.5-10.5:', lucas.ph[~lucas.ph.between(2.5, 10.5)].tolist()[:10])
    lucas['ix'], lucas['iy'] = cell_ids(lucas.la, lucas.lo, args.cell)
    lucas['cell'] = lucas.ix.astype(str) + '_' + lucas.iy.astype(str)
    print(f'LUCAS rows: {n_all} total, {n_ro} Romania, {len(lucas)} with valid pH and coordinates ({bad_ph} without pH)')

    # ---------------- cells
    ce = echo.groupby('cell').agg(ph_e=('soil_ph_field', 'median'), n_e=('soil_ph_field', 'size'),
                                  n_coll=('collected_by', 'nunique'), ix=('ix', 'first'), iy=('iy', 'first')).reset_index()
    cl = lucas.groupby('cell').agg(ph_l=('ph', 'median'), n_l=('ph', 'size')).reset_index()
    m = ce.merge(cl, on='cell', how='inner')
    m['cx'] = (m.ix + 0.5) * args.cell; m['cy'] = (m.iy + 0.5) * args.cell      # cell centres relative to the grid origin
    m['blk'] = np.floor(m.cx / args.block).astype(int).astype(str) + '_' + np.floor(m.cy / args.block).astype(int).astype(str)

    desc = dict(lucas_rows_file=int(n_all), lucas_rows_RO=int(n_ro), lucas_without_valid_ph=int(bad_ph),
                lucas_depth_counts={str(k): int(v) for k, v in depth_counts.items()}, lucas_coord_columns=[lat, lon],
                echo_samples=int(len(echo)), echo_cells=int(len(ce)), lucas_points_RO=int(len(lucas)),
                lucas_cells=int(len(cl)), lucas_density_per_100km2=float(len(lucas) / ROMANIA_AREA_KM2 * 100),
                lucas_km2_per_point=float(ROMANIA_AREA_KM2 / len(lucas)), cells_with_both=int(len(m)),
                lucas_median_all=float(lucas.ph.median()), echo_median_all=float(echo.soil_ph_field.median()))
    sets = {'main': m[(m.n_e >= args.min_echo) & (m.n_l >= args.min_lucas)],
            'S1_ge2_collectors': m[(m.n_e >= args.min_echo) & (m.n_l >= args.min_lucas) & (m.n_coll >= 2)],
            'S2_ge3echo_ge2lucas': m[(m.n_e >= 3) & (m.n_l >= 2)]}
    res = {'descriptive': desc, 'cell_km': args.cell / 1000, 'block_km': args.block / 1000, 'classes': CLASS_LABELS, 'boot': args.boot, 'seed': args.seed}
    for k, d in sets.items():
        if len(d) < 8:
            print(f'WARNING: set {k} has only {len(d)} cells; skipped'); continue
        res[k] = analyse(d, args.boot, rng)
    if 'main' not in res: die('fewer than 8 matched cells in the main set; check coordinates/CRS and the cell size.')
    sets['main'].drop(columns=['ix', 'iy']).to_csv(os.path.join(args.out, 'lucas_cells_matched.csv'), index=False)
    json.dump(res, open(os.path.join(args.out, 'lucas_summary.json'), 'w'), indent=2)

    rows = []
    for k in [s for s in sets if s in res]:
        r = res[k]
        for lab, key, d in [('Spearman rho (cell medians)', 'rho', 2), ('Weighted kappa, 3 pH classes', 'kappa', 2),
                            ('Same class, %', 'same', 0), ('Same or adjacent class, %', 'within1', 0),
                            ('Median ECHO - LUCAS (pH units)', 'meddiff', 2), (f'Spearman rho within {args.block/1000:.0f}-km blocks', 'rho_within', 2)]:
            rows.append(dict(set=k, n_cells=r['n_cells'], metric=lab, estimate=r[key]['est'], ci_low=r[key]['lo'], ci_high=r[key]['hi']))
    pd.DataFrame(rows).to_csv(os.path.join(args.out, 'lucas_table.csv'), index=False)
    make_figure(sets['main'], res, args.out, args.simulated)

    # ---------------- ready-to-paste text
    r = res['main']; d = desc
    txt = []
    txt.append('METHODS (to adapt): LUCAS Topsoil 2018 points for Romania with a valid pH(H2O) value (n = %d; one per %.0f km2; coordinates: %s) and the ECHO samples were assigned to 30 x 30 km cells of the %s grid. Cells with at least %d ECHO samples and %d LUCAS point were retained. Agreement was expressed as the Spearman correlation between the ECHO and LUCAS cell medians, as the linearly weighted kappa for three pH classes (<6.0, 6.0-7.5, >7.5) fixed beforehand, and as the share of cells in the same class; intervals are 95%% percentile intervals from a bootstrap that resamples 100-km blocks. To separate local from regional agreement, the correlation was recomputed after centring both medians within 100-km blocks.'
               % (d['lucas_points_RO'], d['lucas_km2_per_point'], '/'.join(d['lucas_coord_columns']), ('ETRS89-LAEA' if args.crs == 3035 else 'Stereo 70' if args.crs == 3844 else 'EPSG:%d' % args.crs), args.min_echo, args.min_lucas))
    txt.append('RESULTS (to adapt): %d cells contained both ECHO and LUCAS data (%d ECHO samples, %d LUCAS points; %d of %d ECHO cells). Cell medians of ECHO and LUCAS pH were correlated (Spearman rho = %s; n = %d cells). Weighted kappa for the three pH classes was %s; %s%% of the cells fell in the same class and %s%% in the same or an adjacent class. The median difference between ECHO and LUCAS cell medians was %s pH units (descriptive, not an error estimate). After centring within 100-km blocks the correlation was %s.'
               % (r['n_cells'], r['n_echo_samples'], r['n_lucas_points'], r['n_cells'], d['echo_cells'], ci(r['rho']), r['n_cells'],
                  ci(r['kappa']), fmt(r['same']['est'], 0), fmt(r['within1']['est'], 0), ci(r['meddiff']), ci(r['rho_within'])))
    for k in ('S1_ge2_collectors', 'S2_ge3echo_ge2lucas'):
        if k in res:
            s = res[k]
            txt.append(f"SENSITIVITY {k}: n = {s['n_cells']} cells; rho = {ci(s['rho'])}; kappa = {ci(s['kappa'])}.")
    open(os.path.join(args.out, 'lucas_text.txt'), 'w', encoding='utf-8').write('\n\n'.join(txt) + '\n')
    print('\n'.join(txt)); print(f'\nOutputs written to {args.out}/')

if __name__ == '__main__':
    main()
