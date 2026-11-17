#!/usr/bin/env python3
"""Summary statistics behind the tables in README.md.

    python scripts/summary_statistics.py --workdir ./work [--season 2025-26] [--out results]

Reads work/storms/ndfd_<OFFICE>.pkl (made by the `storms` stage of nws_snowfall_verification.py) and writes CSV files:
    office_summary.csv            area-wide metrics per office, full record and one season (with storm-bootstrap 90% interval for mean error)
    error_by_observed_amount.csv  mean error / hit / high / low rates by observed storm-total bin   (conditioning on O: see README, section 5.4)
    error_by_forecast_amount.csv  the same by forecast amount bin                                    (conditioning on F: what a forecast user needs)
    error_by_amount_pooled.csv    both amount tables pooled over all six offices (case-weighted)
    error_by_terrain.csv          mean error and percent bias by coast / lowland / hills / mountains
    largest_storms.csv            the five largest observed storms per office with storm-mean F and O
All metrics use exactly the masks of the figures (office CWA cells; radar beam-blockage sectors excluded).
"""
import argparse, csv, datetime as dt, importlib.util, os, pathlib, pickle, sys
import numpy as np
from scipy.ndimage import distance_transform_edt, map_coordinates

ROOT = pathlib.Path(__file__).resolve().parents[1]
ap = argparse.ArgumentParser(); ap.add_argument("--workdir", default="work"); ap.add_argument("--season", default="2025-26"); ap.add_argument("--out", default=str(ROOT / "results"))
args = ap.parse_args()
os.environ["NWSNOW_WORKDIR"] = str(pathlib.Path(args.workdir).resolve())
os.environ.setdefault("NWSNOW_CWA_FILE", str(ROOT / "data" / "cwa_new_england.geojson"))
spec = importlib.util.spec_from_file_location("nsv", ROOT / "nws_snowfall_verification.py"); nsv = importlib.util.module_from_spec(spec); spec.loader.exec_module(nsv)
config = nsv.config; OUT = pathlib.Path(args.out); OUT.mkdir(parents=True, exist_ok=True)
rng = np.random.default_rng(2026)

def terrain_vars(CL, CA):
    z = np.load(config.DEM); lon, lat, elev = z["lon"], z["lat"], z["elev"]
    if lat[1] < lat[0]: lat = lat[::-1]; elev = elev[::-1]
    ix = (CL[None, :] - lon[0]) / (lon[1] - lon[0]) + 0 * CA[:, None]; iy = (CA[:, None] - lat[0]) / (lat[1] - lat[0]) + 0 * CL[None, :]
    EL = np.maximum(map_coordinates(elev, [iy, ix], order=1, mode="nearest"), 0)
    sea = elev <= 0.5; la0 = np.searchsorted(lat, max(CA.min() - 3, lat[0])); la1 = np.searchsorted(lat, min(CA.max() + 3, lat[-1]))
    lo0 = np.searchsorted(lon, CL.min() - 3); lo1 = np.searchsorted(lon, min(CL.max() + 3, lon[-1]))
    sub = (slice(la0, la1), slice(lo0, lo1)); ds = distance_transform_edt(~sea[sub]) * (abs(lat[1] - lat[0]) * 111.0)
    return EL, map_coordinates(ds, [iy - la0, ix - lo0], order=1, mode="nearest")

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
    keys = sorted(R)
    F = np.stack([R[k]["fc"] for k in keys]); T = np.stack([R[k]["nsa"] for k in keys])
    return dict(R=R, keys=keys, F=F, T=T, E=F - T, CM=CM, XM=XM, CL=CL, CA=CA, ok=np.isfinite(F - T) & CM[None])

def tol_of(T): return np.maximum(2.0, 0.30 * T)

def office_row(off, per, d):
    F, T, E, ok, CM, R, keys = d["F"], d["T"], d["E"], d["ok"], d["CM"], d["R"], d["keys"]; NS = len(keys)
    hit = np.abs(E) <= tol_of(T); over = E > tol_of(T); under = E < -tol_of(T); heavy = ok & (np.maximum(F, T) >= 4.0)
    fo = np.array([F[i][ok[i]].mean() for i in range(NS)]); ob = np.array([T[i][ok[i]].mean() for i in range(NS)]); sm = fo - ob
    boot = np.array([sm[rng.integers(0, NS, NS)].mean() for _ in range(5000)]); lo, hi = np.percentile(boot, [5, 95])
    with np.errstate(all="ignore"): mE = np.nanmean(np.where(ok, E, np.nan), 0)
    v = mE[np.isfinite(mE) & CM]
    return dict(office=off, period=per, storms=NS, first_storm=min(R[k]["w0"] for k in keys).strftime("%Y-%m-%d"), last_storm=max(R[k]["w1"] for k in keys).strftime("%Y-%m-%d"),
                storms_12h_lead=sum(R[k]["lead"] == 12 for k in keys), storms_6h_lead=sum(R[k]["lead"] == 6 for k in keys), cells_used=int(CM.sum()), cells_excluded_radar=int(d["XM"].sum()),
                mean_forecast_in=F[ok].mean(), mean_observed_in=T[ok].mean(), mean_error_in=E[ok].mean(), mean_error_ci90_lo=lo, mean_error_ci90_hi=hi,
                mae_in=np.abs(E[ok]).mean(), rmse_in=np.sqrt((E[ok] ** 2).mean()), percent_bias=100 * E[ok].sum() / T[ok].sum(), normalized_mae=100 * np.abs(E[ok]).sum() / T[ok].sum(),
                hit_rate=100 * hit[ok].mean(), high_bias_rate=100 * over[ok].mean(), low_bias_rate=100 * under[ok].mean(),
                heavy_cases=int(heavy.sum()), heavy_hit_rate=100 * hit[heavy].mean(), heavy_high_rate=100 * over[heavy].mean(), heavy_low_rate=100 * under[heavy].mean(),
                storm_mean_corr=np.corrcoef(fo, ob)[0, 1] if NS > 2 else float("nan"), storms_within_band_pct=100 * np.mean(np.abs(sm) <= np.maximum(0.25 * ob, 1.0)),
                storms_high_gt05=int((sm > 0.5).sum()), storms_low_lt05=int((sm < -0.5).sum()),
                area_pct_mean_err_gt_p05=100 * np.mean(v > 0.5), area_pct_mean_err_lt_m05=100 * np.mean(v < -0.5))

BINS = [0, 1, 2, 4, 6, 9, 12, 1000]; LAB = ["<1", "1-2", "2-4", "4-6", "6-9", "9-12", "12+"]
def amount_rows(off, per, d, by):
    X = (d["T"] if by == "O" else d["F"])[d["ok"]]; E = d["E"][d["ok"]]; T = d["T"][d["ok"]]; tol = tol_of(T); rows = []
    for i, l in enumerate(LAB):
        s = (X >= BINS[i]) & (X < BINS[i + 1])
        if s.sum() >= 50: rows.append(dict(office=off, period=per, bin_in=l, cases=int(s.sum()), mean_error_in=E[s].mean(), mean_forecast_in=d["F"][d["ok"]][s].mean(), mean_observed_in=T[s].mean(),
                                          hit_rate=100 * (np.abs(E[s]) <= tol[s]).mean(), high_bias_rate=100 * (E[s] > tol[s]).mean(), low_bias_rate=100 * (E[s] < -tol[s]).mean()))
    return rows

def terrain_rows(off, per, d):
    EL, DK = terrain_vars(d["CL"], d["CA"]); CM = d["CM"]; rows = []
    for n, m in {"coast (<15 km from ocean)": CM & (DK < 15), "inland lowland (<200 m)": CM & (DK >= 15) & (EL < 200),
                 "hills (200-600 m)": CM & (DK >= 15) & (EL >= 200) & (EL < 600), "mountains (>=600 m)": CM & (DK >= 15) & (EL >= 600)}.items():
        if m.sum() >= 15:
            mm = d["ok"] & m[None]; rows.append(dict(office=off, period=per, terrain_class=n, cells=int(m.sum()), mean_observed_in=d["T"][mm].mean(), mean_error_in=d["E"][mm].mean(), percent_bias=100 * d["E"][mm].sum() / d["T"][mm].sum()))
    return rows

def storm_rows(off, d):
    ok, F, T = d["ok"], d["F"], d["T"]; fo = np.array([F[i][ok[i]].mean() for i in range(len(d["keys"]))]); ob = np.array([T[i][ok[i]].mean() for i in range(len(d["keys"]))])
    return [dict(office=off, storm_start=d["R"][d["keys"][i]]["label"], window=d["R"][d["keys"][i]]["window"], observed_storm_mean_in=ob[i], forecast_storm_mean_in=fo[i], error_in=fo[i] - ob[i], forecast_run=d["R"][d["keys"][i]]["issuance"]) for i in np.argsort(-ob)[:5]]

def write(name, rows):
    if not rows: return
    with open(OUT / name, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0])); w.writeheader()
        for r in rows: w.writerow({k: (f"{v:.3f}" if isinstance(v, (float, np.floating)) else v) for k, v in r.items()})
    print("wrote", OUT / name, len(rows), "rows")

summ, byO, byF, ter, big = [], [], [], [], []; pool = {"full": {"O": [], "F": []}, args.season: {"O": [], "F": []}}
for off in config.ORDER:
    for per, season in (("full record", None), (f"{args.season} season", args.season)):
        d = load(off, season)
        if d is None: continue
        summ.append(office_row(off, per, d)); byO += amount_rows(off, per, d, "O"); byF += amount_rows(off, per, d, "F"); ter += terrain_rows(off, per, d)
        if season is None: big += storm_rows(off, d)
def pooled(rows, by):
    out = []
    for per in sorted({r["period"] for r in rows}):
        for l in LAB:
            rs = [r for r in rows if r["period"] == per and r["bin_in"] == l]
            if not rs: continue
            n = sum(r["cases"] for r in rs); w = lambda k: sum(r["cases"] * r[k] for r in rs) / n
            out.append(dict(binned_by=by, period=per, bin_in=l, cases=n, mean_error_in=w("mean_error_in"), mean_forecast_in=w("mean_forecast_in"), mean_observed_in=w("mean_observed_in"),
                            hit_rate=w("hit_rate"), high_bias_rate=w("high_bias_rate"), low_bias_rate=w("low_bias_rate")))
    return out
write("error_by_amount_pooled.csv", pooled(byO, "observed") + pooled(byF, "forecast"))
write("office_summary.csv", summ); write("error_by_observed_amount.csv", byO); write("error_by_forecast_amount.csv", byF); write("error_by_terrain.csv", ter); write("largest_storms.csv", big)
