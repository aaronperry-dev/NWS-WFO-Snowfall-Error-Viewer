#!/usr/bin/env python3
"""Export the data the web viewer needs: one raw float32 .bin grid per office (mean snowfall error), a manifest, and CWA outlines.

    python scripts/export_viewer_data.py --workdir ./work [--season 2025-26] [--out viewer/data]

--season defaults to the most recent season present in the storm files (a season runs Jul 1 .. Jun 30, e.g. 2025-26).
Reads work/storms/ndfd_<OFFICE>.pkl (made by the `storms` stage of nws_snowfall_verification.py) and writes, under --out:

    manifest.json                       what exists, grid geometry per file, per-office statistics (merged across runs, so seasons accumulate)
    <season>/<OFFICE>_mean_error.bin    float32 little-endian, row-major, row 0 = NORTH edge, no header (geometry is in the manifest)
    cwa.geojson                         simplified CWA outlines for the six offices (EPSG:4326)

Grid values, exactly as in the printed maps (figure stage e): mean over the season's storms of (forecast - observed), inches, for cells
whose center lies inside the CWA and that have valid data in at least 90 % of the storms.
    NaN     outside the CWA / too few valid storms
    -999    radar beam-blockage sector (config.RADAR_EXCLUSIONS): excluded from every statistic
"""
import argparse, csv, datetime as dt, importlib.util, json, os, pathlib, pickle, sys
import numpy as np, shapely

ROOT = pathlib.Path(__file__).resolve().parents[1]
ap = argparse.ArgumentParser(); ap.add_argument("--workdir", default="work"); ap.add_argument("--season", default=None); ap.add_argument("--out", default=str(ROOT / "viewer" / "data"))
args = ap.parse_args()
os.environ["NWSNOW_WORKDIR"] = str(pathlib.Path(args.workdir).resolve())
os.environ.setdefault("NWSNOW_CWA_FILE", str(ROOT / "data" / "cwa_new_england.geojson"))
spec = importlib.util.spec_from_file_location("nsv", ROOT / "nws_snowfall_verification.py"); nsv = importlib.util.module_from_spec(spec); spec.loader.exec_module(nsv)
config = nsv.config; OUT = pathlib.Path(args.out); EXCLUDED = -999.0; MIN_VALID = 0.9; SMALL_SAMPLE = 10
rng = np.random.default_rng(2026)

def season_of(t): return t.year if t.month >= 7 else t.year - 1
def label(season): y = int(season[:4]); return f"{y}–{str(y + 1)[2:]}"

def load(off, season=None):
    D = pickle.load(open(config.STORMS / f"ndfd_{off}.pkl", "rb")); R = D["R"]; CL, CA, CM = D["CLON"], D["CLAT"], D["CMASK"].copy()
    if season:
        y = int(season[:4]); R = {k: v for k, v in R.items() if dt.datetime(y, 7, 1) <= v["w0"] < dt.datetime(y + 1, 7, 1)}
    if not R: return None
    XM = np.zeros(CM.shape, bool)
    if off in config.RADAR_EXCLUSIONS:
        nm, la0, lo0, a0, a1, r0 = config.RADAR_EXCLUSIONS[off]; LON, LAT = np.meshgrid(CL, CA)
        dy = (LAT - la0) * 111.0; dx = (LON - lo0) * 111.0 * np.cos(np.radians(la0)); az = (np.degrees(np.arctan2(dx, dy)) + 360) % 360
        XM = CM & (az >= a0) & (az <= a1) & (np.hypot(dx, dy) > r0); CM = CM & ~XM
    keys = sorted(R); F = np.stack([R[k]["fc"] for k in keys]); T = np.stack([R[k]["nsa"] for k in keys]); E = F - T
    return dict(R=R, keys=keys, F=F, T=T, E=E, CM=CM, XM=XM, CL=CL, CA=CA, ok=np.isfinite(E) & CM[None])

def grid_and_stats(d):
    E, F, T, ok, CM, XM, keys, R = d["E"], d["F"], d["T"], d["ok"], d["CM"], d["XM"], d["keys"], d["R"]; NS = len(keys)
    n = ok.sum(0).astype(float); n[n < MIN_VALID * NS] = np.nan
    with np.errstate(all="ignore"): g = np.nanmean(np.where(ok, E, np.nan), 0)
    g[np.isnan(n)] = np.nan; g[~CM] = np.nan; g = g.astype("<f4"); g[XM] = EXCLUDED
    tol = np.maximum(2.0, 0.30 * T); sm = np.array([E[i][ok[i]].mean() for i in range(NS)])
    boot = np.array([sm[rng.integers(0, NS, NS)].mean() for _ in range(5000)]); lo, hi = np.percentile(boot, [5, 95])
    v = g[np.isfinite(g) & (g != EXCLUDED)]
    st = dict(storms=NS, first_storm=min(R[k]["w0"] for k in keys).strftime("%Y-%m-%d"), last_storm=max(R[k]["w1"] for k in keys).strftime("%Y-%m-%d"),
              mean_error=float(E[ok].mean()), mean_error_ci90=[float(lo), float(hi)], mae=float(np.abs(E[ok]).mean()), rmse=float(np.sqrt((E[ok] ** 2).mean())),
              percent_bias=float(100 * E[ok].sum() / T[ok].sum()), hit_rate=float(100 * (np.abs(E[ok]) <= tol[ok]).mean()),
              mean_forecast=float(F[ok].mean()), mean_observed=float(T[ok].mean()),
              area_pct_gt_p05=float(100 * np.mean(v > 0.5)), area_pct_lt_m05=float(100 * np.mean(v < -0.5)),
              cells=int(CM.sum()), cells_excluded_radar=int(XM.sum()), small_sample=bool(NS < SMALL_SAMPLE))
    return g, st

def ci_from_results(off, season, st):
    """Reuse the 90 % storm-bootstrap interval already published in results/office_summary.csv (written by scripts/summary_statistics.py) so the viewer
    and the README quote identical numbers.  Used only if that row describes the same storms and mean error; otherwise the interval computed here stays."""
    f = ROOT / "results" / "office_summary.csv"
    if not f.exists(): return None
    for r in csv.DictReader(open(f)):
        if r["office"] == off and r["period"] == f"{season} season" and int(r["storms"]) == st["storms"] and abs(float(r["mean_error_in"]) - st["mean_error"]) < 1e-3:
            return [float(r["mean_error_ci90_lo"]), float(r["mean_error_ci90_hi"])]
    return None

# ---- which season
if args.season: season = args.season
else:
    last = max(max(v["w0"] for v in pickle.load(open(config.STORMS / f"ndfd_{o}.pkl", "rb"))["R"].values()) for o in config.ORDER if (config.STORMS / f"ndfd_{o}.pkl").exists())
    season = f"{season_of(last)}-{str(season_of(last) + 1)[2:]}"
print("season", season)
(OUT / season).mkdir(parents=True, exist_ok=True)

# ---- grids + stats
offices = {}
for off in config.ORDER:
    if not (config.STORMS / f"ndfd_{off}.pkl").exists(): print(off, ": no storm file, skipped"); continue
    d = load(off, season)
    if d is None: print(off, f": no storms in {season}, skipped"); continue
    g, st = grid_and_stats(d); st["mean_error_ci90"] = ci_from_results(off, season, st) or st["mean_error_ci90"]; full = load(off); st["full_record_mean_error"] = float(full["E"][full["ok"]].mean()); st["full_record_storms"] = len(full["keys"])
    ny, nx = g.shape; dx = abs(float(d["CL"][1] - d["CL"][0])); dy = abs(float(d["CA"][1] - d["CA"][0]))
    assert d["CA"][0] > d["CA"][-1], "row 0 must be the northern edge"
    fn = f"{season}/{off}_mean_error.bin"; g.tofile(OUT / fn); assert (OUT / fn).stat().st_size == nx * ny * 4
    offices[off] = dict(file=fn, nx=nx, ny=ny, west=float(d["CL"][0] - dx / 2), north=float(d["CA"][0] + dy / 2), dx=dx, dy=dy, stats=st)
    print(f"{off}: {ny}x{nx} cells, {st['storms']} storms, mean error {st['mean_error']:+.2f} in -> {fn} ({(OUT / fn).stat().st_size / 1024:.0f} KB)")

# ---- manifest (merged with what is already there so seasons accumulate)
mf = OUT / "manifest.json"; M = json.load(open(mf)) if mf.exists() else {}
M.update(version=1, units="in", byte_order="little-endian float32, row-major, row 0 = north edge, no header", nodata="NaN", excluded=EXCLUDED,
         layer=dict(id="mean_error", label="Mean error (forecast − observed)", units="inches of snow per storm", vmin=-2.0, vmax=2.0),
         offices={o: dict(name=config.OFFICES[o]["name"], color=config.OFFICES[o]["color"], point=list(config.OFFICES[o]["pt"])) for o in config.ORDER},
         radar_exclusions={o: dict(radar=v[0], azimuth=[v[3], v[4]], min_range_km=v[5]) for o, v in config.RADAR_EXCLUSIONS.items()},
         min_valid_storm_fraction=MIN_VALID, small_sample_storms=SMALL_SAMPLE)
M.setdefault("seasons", {})[season] = dict(label=f"{label(season)} season", window=f"Jul 1 {season[:4]} – Jun 30 {int(season[:4]) + 1}", offices=offices)
M["default_season"] = max(M["seasons"])
json.dump(M, open(mf, "w"), indent=1); print("wrote", mf)

# ---- CWA outlines (simplified like the printed maps: islets dropped, 0.008 deg tolerance)
feats = []
for off, g in nsv.load_cwa().items():
    if off not in config.OFFICES: continue
    o = shapely.set_precision(nsv.outline_geom(g), 0.0001); c = g.representative_point()
    feats.append(dict(type="Feature", properties=dict(cwa=off, name=config.OFFICES[off]["name"], label=[round(c.x, 4), round(c.y, 4)]), geometry=json.loads(shapely.to_geojson(o))))
json.dump(dict(type="FeatureCollection", features=feats), open(OUT / "cwa.geojson", "w"), separators=(",", ":")); print("wrote", OUT / "cwa.geojson", (OUT / "cwa.geojson").stat().st_size // 1024, "KB")
