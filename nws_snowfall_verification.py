#!/usr/bin/env python3
"""
NWS New England snowfall-forecast verification: one self-contained script, start to finish.

Verifies the official WFO-edited NDFD snowfall forecasts (NCEI archive, June 2018 on) against NOHRSC sfav2 6-h snowfall analyses for the
New England forecast offices (CAR, GYX, BTV, ALY, BOX, OKX) and draws native-resolution (unsmoothed, ~4 km) maps over terrain shading.

    pip install numpy scipy matplotlib cartopy shapely pyproj pygrib tifffile pillow requests     # pygrib needs the ECCODES library

    python nws_snowfall_verification.py --workdir ./work                          # everything, all six offices
    python nws_snowfall_verification.py --workdir ./work --offices GYX BTV        # a subset of offices
    python nws_snowfall_verification.py --workdir ./work --stages figures --force # just redraw the figures
    python nws_snowfall_verification.py --workdir ./work --stages figures --season 2025-26   # one winter only -> figures/season_2025-26/
    python nws_snowfall_verification.py --workdir ./work --start 2022-11-01 --end 2023-04-30 --offices GYX BTV   # quick test, ~3 min
    python nws_snowfall_verification.py --workdir ./work --jobs 4 --nohrsc-rate 4 --ncei-rate 2                   # gentler on the servers
    python nws_snowfall_verification.py radar-test GYX -70.257 43.891 [out.png]    # radar beam-blockage check (see RADAR_EXCLUSIONS)

Stages (run in this order; each is skipped when its outputs already exist unless --force; every stage is resume-safe):
    terrain  AWS terrain tiles                                  -> work/terrain/dem.npz
    cube     NOHRSC 6-h snowfall, New England window            -> work/nohrsc_cube/*.npy        (~1.7 GB for the full record; the slow stage)
    storms   storm detection + NDFD forecasts regridded         -> work/storms/ndfd_<OFFICE>.pkl (downloads ~2 GB of NDFD GRIB2 for the full record)
    figures  PNG maps and charts                                -> work/figures/                 (~30-50 s per office; offices run in parallel)

Figures per office (--figures e,a,m,b,g,c,d):
    e  _0_mean_error_map           large single-panel Mean Error (F - O, in) averaged over storms
    a  _1_bias_MAE_maps            Mean Error, MAE, Percent Bias, Normalized MAE
    m  _6_mean_F_O_error_maps      mean forecast, mean observed, mean error
    b  _2_hit_rate_maps            Hit Rate, Forecast Tendency (high / low bias), Reliability Category
    g  _5_largest_storms_F_O_error the four largest storms: forecast, observed, error
    c  _3_error_charts             hit / high / low rates by snowfall amount and terrain; F vs O
    d  _4_bias_terrain_coast       Percent Bias vs elevation and distance to ocean
plus regional_mean_error_map.png (all offices except REGIONAL_EXCLUDE on one map and one color scale).  With --season the default is e,a,m.

How the numbers are made
    Forecast  official WFO-edited NDFD "Snowfall Accumulation" (WMO YSUZ98, KWBN, CONUS 2.5 km, 6-h periods), NCEI THREDDS archive.  Not NBM.
              Forecast = latest run issued >= 12 h (fallback 6 h) before the storm window opens; regridded to NOHRSC cells by 5x5 sub-sample
              bilinear averaging (one sparse-matrix multiply per message).
    Observed  NOHRSC sfav2 6-h snowfall (0.04 deg).
    Storms    CWA-mean 6-h series >= 0.08 in, gaps <= 12 h merged, cluster total >= 1.5 in; events longer than 48 h split at the quietest 6 h.
    Metrics   error E = F - O.  Percent Bias = 100*sum(E)/sum(O); Normalized MAE = 100*sum|E|/sum(O) (cells with mean O < 1 in masked).
              Hit = |F - O| <= max(2 in, 30% of O); otherwise high bias (F > O) or low bias (F < O).  Heavy-snow cases: max(F, O) >= 4 in,
              at least 8 per cell.  Reliability: reliable if hit >= 60%; high / low bias if hit < 60% and one tendency exceeds the other by
              >= 20 points; otherwise mixed.

Cartography (every map)
    Lambert conformal conic (central meridian 71.5W, standard parallels 42N / 46N).  Layer order: terrain hillshade, water, data (native
    cells, masked to the CWA), then boundaries above the data: international boundary (solid black, white casing), state / province lines
    (dashed gray), coastline (thin blue; islets < ~5 km2 generalized away), focus CWA (office color, white halo), neighboring CWAs (dash-dot,
    labelled).  Graticule labelled inside the neatline, state / water-body names outside the data, reference cities, collision-free labels,
    a legend limited to the symbols present, dual mi / km scale bar, north arrow rotated for grid convergence, data credit.  One fixed
    diverging color scale per quantity for every office (+/-2 in for mean error); colorbars get arrow ends when values are clipped.
    Radar beam-blockage sectors are hatched, excluded from every statistic and listed in the legend.

Network politeness: keep-alive connections, per-host connection and requests-per-second caps (NET_LIMITS), exponential backoff, Retry-After,
automatic slow-down on 429/503, a User-Agent, request / MB / retry counts printed at the end of each stage.  Everything downloaded is cached.

Inputs: the CWA polygons are read from --cwa-file, else data/cwa_new_england.geojson next to this script, else <workdir>/cwa_new_england.geojson;
if none exists they are downloaded once from the NWS reference map service (mapservices.weather.noaa.gov) and cached in the work directory.

Environment overrides: NWSNOW_WORKDIR, NWSNOW_START, NWSNOW_END, NWSNOW_NOHRSC_RATE, NWSNOW_NCEI_RATE, NWSNOW_USER_AGENT, NWSNOW_CWA_FILE.

Caveats: NOHRSC is an analysis, not ground truth.  NDFD archive coverage begins June 2018.  Disk: ~6-8 GB for the full record.  RAM: ~2.5-3 GB per
figure process.  The first figure run downloads Natural Earth 10 m shapefiles (cartopy, one time).  Forecaster staffing changes are not accounted for.
To add an office: add it to OFFICES / ORDER below (its CWA polygon is then fetched automatically).
"""
import argparse, glob, io, json, math, os, pathlib, pickle, random, re, shutil, subprocess, sys, textwrap, threading, time, types
import datetime as dt
from collections import Counter
from concurrent.futures import ThreadPoolExecutor, ProcessPoolExecutor
from urllib.parse import urlparse, urlencode

try:
    import numpy as np
    import scipy.sparse as sp
    from scipy.ndimage import gaussian_filter, zoom, distance_transform_edt, map_coordinates
    from scipy.spatial.distance import cdist
    import requests
    from requests.adapters import HTTPAdapter
    from PIL import Image
    import shapely
    from shapely import contains_xy, make_valid
    from shapely.geometry import shape
    from pyproj import Transformer
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    import matplotlib.patheffects as pe
    from matplotlib.colors import LightSource, LinearSegmentedColormap, ListedColormap, BoundaryNorm
    from matplotlib.patches import Patch, Rectangle
    from matplotlib.lines import Line2D
    from matplotlib.collections import LineCollection
    from matplotlib.figure import Figure
    import cartopy.crs as ccrs
    import cartopy.feature as cfeature
    import cartopy.io.shapereader as shp
    from cartopy.mpl.geoaxes import GeoAxes
except ImportError as _e:
    sys.exit(f"Missing Python package ({_e}).\nInstall with:  pip install numpy scipy matplotlib cartopy shapely pyproj pygrib tifffile pillow requests")

SELF = pathlib.Path(__file__).resolve()

# ======================================================================================================================================
# SETTINGS  (edit here)
# ======================================================================================================================================
# ---- offices (name, map color, office location)
OFFICES = {
    "BOX": dict(name="Boston/Norton",  color="#E69F00", pt=(-71.18, 41.97)),
    "GYX": dict(name="Gray/Portland",  color="#0072B2", pt=(-70.26, 43.89)),
    "CAR": dict(name="Caribou",        color="#009E73", pt=(-68.02, 46.87)),
    "BTV": dict(name="Burlington",     color="#CC79A7", pt=(-73.15, 44.47)),
    "ALY": dict(name="Albany",         color="#D55E00", pt=(-73.75, 42.65)),
    "OKX": dict(name="New York/Upton", color="#56B4E9", pt=(-72.86, 40.87)),
}
ORDER = ["CAR", "GYX", "BTV", "ALY", "BOX", "OKX"]
REGIONAL_EXCLUDE = ("OKX",)          # offices left off the combined regional map (small sample: only a few dozen storms in the record)

# ---- study period (cold-season months only); NDFD archive begins June 2018.  Override with --start / --end.
DEFAULT_START, DEFAULT_END = "2018-11-01", "2026-04-30"
WINTER_MONTHS = (11, 12, 1, 2, 3, 4)

# ---- NOHRSC sfav2 CONUS 6-h snowfall grid: 1500 x 850 cells of 0.04 deg, upper-left corner (-126, 55)
NOHRSC_X0, NOHRSC_Y0, NOHRSC_STEP = -126.0, 55.0, 0.04
CUBE_ROWS = (175, 388)               # rows  175..387  (lat 48.0 .. 39.5)
CUBE_COLS = (1150, 1500)             # cols 1150..1499 (lon -80.0 .. -66.0)   -> covers all six CWAs

# ---- storm detection (CWA-mean 6-h NOHRSC snowfall)
THR_SERIES = 0.08                    # in per 6 h: a period counts as "snowing" at or above this CWA-mean value
THR_TOTAL  = 1.5                     # in: minimum CWA-mean storm total to keep a storm
GAP_HOURS  = 12                      # periods closer than this are merged into one storm
MAX_HOURS  = 48                      # longer events are split at their quietest 6-h period
LEADS      = (12, 6)                 # forecast issued at least this many hours before the window opens (12 h, then 6 h fallback)

# ---- radar beam-blockage exclusions: cells of the NOHRSC analysis that are artifacts of radar geometry
# office -> (radar name, radar lat, radar lon, azimuth start deg, azimuth end deg, minimum range km)
# Found with the `radar-test` command; edit or add entries and re-run the figures stage.
RADAR_EXCLUSIONS = {
    "BTV": ("KCXX Burlington", 44.511, -73.166, 30.0, 160.0, 15.0),
    "GYX": ("KGYX Gray",       43.891, -70.257, 69.0,  82.0, 15.0),
}

# ---- network politeness.  None of these hosts publishes a numeric rate limit, so the defaults are conservative.  A 429/503 halves the
# pacing automatically (and Retry-After is honored), then slowly recovers.  Rates can also be set with --nohrsc-rate / --ncei-rate.
NET_LIMITS = {   # host: (max simultaneous connections, max requests per second; 0 = unlimited)
    "www.nohrsc.noaa.gov": (6, 8.0),
    "www.ncei.noaa.gov":   (4, 3.0),
    "s3.amazonaws.com":    (8, 0),
}
NET_DEFAULT = (4, 5.0)
CWA_URL = "https://mapservices.weather.noaa.gov/static/rest/services/nws_reference_maps/nws_reference_map/MapServer/1/query"

# ======================================================================================================================================
# CONFIGURATION OBJECT  (paths and period depend on the work directory, so they are (re)built from the environment)
# ======================================================================================================================================
config = types.SimpleNamespace()

def _init_config():
    c = config; env = os.environ.get
    for k in ("OFFICES", "ORDER", "REGIONAL_EXCLUDE", "WINTER_MONTHS", "NOHRSC_X0", "NOHRSC_Y0", "NOHRSC_STEP", "CUBE_ROWS", "CUBE_COLS",
              "THR_SERIES", "THR_TOTAL", "GAP_HOURS", "MAX_HOURS", "LEADS", "RADAR_EXCLUSIONS", "NET_DEFAULT"):
        setattr(c, k, globals()[k])
    c.WORK = pathlib.Path(env("NWSNOW_WORKDIR", "work")).resolve()
    c.TERRAIN = c.WORK / "terrain"; c.DEM = c.TERRAIN / "dem.npz"           # AWS terrain tiles + dem.npz
    c.CUBE = c.WORK / "nohrsc_cube"                                        # NOHRSC 6-h snowfall, New England window, one .npy per period
    c.NDFD = c.WORK / "ndfd"; c.CATALOG = c.NDFD / "catalog"; c.RUNS = c.NDFD / "runs"   # THREDDS listings, downloaded GRIB2 run files
    c.STORMS = c.WORK / "storms"                                           # ndfd_<OFFICE>.pkl : storm list with forecast and observed grids
    c.FIGS = c.WORK / "figures"; c.LOGS = c.WORK / "logs"
    c.START = dt.datetime.strptime(env("NWSNOW_START", DEFAULT_START), "%Y-%m-%d")
    c.END = dt.datetime.strptime(env("NWSNOW_END", DEFAULT_END), "%Y-%m-%d")
    c.USER_AGENT = env("NWSNOW_USER_AGENT", "nws-snowfall-verification/1.0 (research use; python-requests)")
    c.NET = dict(NET_LIMITS)
    if env("NWSNOW_NOHRSC_RATE"): c.NET["www.nohrsc.noaa.gov"] = (NET_LIMITS["www.nohrsc.noaa.gov"][0], float(env("NWSNOW_NOHRSC_RATE")))
    if env("NWSNOW_NCEI_RATE"): c.NET["www.ncei.noaa.gov"] = (NET_LIMITS["www.ncei.noaa.gov"][0], float(env("NWSNOW_NCEI_RATE")))

_init_config()


# ======================================================================================================================================
# POLITE, RESILIENT HTTP  (used by every download stage)
# keep-alive sessions; per-host caps on connections and requests/second (NET_LIMITS); retries with exponential backoff + jitter on 429/5xx/
# timeouts, honoring Retry-After; a 429/503 doubles the spacing between requests to that host and sustained success speeds it back up;
# User-Agent; per-host request / byte / retry counters (http_summary() at the end of a stage)
# ======================================================================================================================================
RETRY_STATUS = {429, 500, 502, 503, 504}
_tl = threading.local()
_hosts = {}; _hl = threading.Lock()

class Host:
    def __init__(self, name):
        conn, rate = config.NET.get(name, config.NET_DEFAULT)
        self.name = name; self.sem = threading.BoundedSemaphore(conn)
        self.base = 1.0 / rate if rate else 0.0; self.interval = self.base
        self.lock = threading.Lock(); self.next = 0.0; self.pause_until = 0.0; self.streak = 0
        self.requests = 0; self.bytes = 0; self.retries = 0; self.throttled = 0; self.t0 = time.monotonic()
    def pace(self):
        with self.lock:
            now = time.monotonic(); t = max(now, self.next, self.pause_until); self.next = t + self.interval
        if t > now: time.sleep(t - now)
    def ok(self, nbytes):
        with self.lock:
            self.requests += 1; self.bytes += nbytes; self.streak += 1
            if self.streak >= 50 and self.interval > self.base: self.interval = max(self.base, self.interval * 0.8); self.streak = 0
    def slow(self, retry_after=None):
        with self.lock:
            self.throttled += 1; self.streak = 0
            self.interval = min(max(self.interval * 2, 0.25), 5.0)
            self.pause_until = max(self.pause_until, time.monotonic() + (retry_after if retry_after else 5.0))

def _host(name):
    with _hl:
        if name not in _hosts: _hosts[name] = Host(name)
        return _hosts[name]

def _session():
    if not hasattr(_tl, "s"):
        s = requests.Session(); s.headers["User-Agent"] = config.USER_AGENT
        a = HTTPAdapter(pool_connections=4, pool_maxsize=8, max_retries=0); s.mount("https://", a); s.mount("http://", a)
        _tl.s = s
    return _tl.s

def http_get(url, headers=None, ok404=False, tries=6, timeout=(15, 120)):
    """GET with pacing and retries. Returns the Response (200/206), or None for a 404 when ok404. Raises on persistent failure."""
    h = _host(urlparse(url).netloc); last = None
    for k in range(tries):
        h.sem.acquire()
        try:
            h.pace()
            r = _session().get(url, headers=headers, timeout=timeout)
            if r.status_code in (200, 206):
                h.ok(len(r.content)); return r
            if r.status_code == 404 and ok404:
                h.ok(0); return None
            ra = r.headers.get("Retry-After"); ra = float(ra) if ra and ra.replace(".", "", 1).isdigit() else None
            last = f"HTTP {r.status_code}"
            if r.status_code not in RETRY_STATUS: raise RuntimeError(f"{last} for {url}")
            if r.status_code in (429, 503): h.slow(min(ra, 120) if ra else None)
        except (requests.ConnectionError, requests.Timeout, requests.exceptions.ChunkedEncodingError) as e:
            last = e.__class__.__name__
        finally:
            h.sem.release()
        h.retries += 1
        time.sleep(min(60, 2 ** k) * (0.5 + random.random()))
    raise RuntimeError(f"giving up on {url}: {last}")

def http_download(url, out, **kw):
    """GET url and write it atomically to `out`. False if the server has no such file."""
    r = http_get(url, ok404=True, **kw)
    if r is None: return False
    tmp = f"{out}.{os.getpid()}.{threading.get_ident()}.part"
    with open(tmp, "wb") as f: f.write(r.content)
    os.replace(tmp, out); return True

def http_summary():
    for h in _hosts.values():
        dt_ = max(time.monotonic() - h.t0, 1e-6)
        print(f"  net {h.name}: {h.requests} requests, {h.bytes/1e6:.0f} MB, {h.requests/dt_:.1f} req/s avg, "
              f"{h.retries} retries, {h.throttled} throttle responses (pacing now {1/h.interval if h.interval else 0:.1f} req/s)", flush=True)

# ======================================================================================================================================
# CWA POLYGONS
# ======================================================================================================================================
_CWA = {}

def cwa_file():
    """path of the CWA polygon file: NWSNOW_CWA_FILE / --cwa-file, else data/ next to this script, else the cached download in the work dir"""
    for f in (os.environ.get("NWSNOW_CWA_FILE"), SELF.parent / "data" / "cwa_new_england.geojson", config.WORK / "cwa_new_england.geojson"):
        if f and pathlib.Path(f).exists(): return pathlib.Path(f)
    return None

def fetch_cwa_file(out):
    """download the CWA polygons of OFFICES from the NWS reference map service (layer 1, 'NWS County Warning Area'), lightly simplified"""
    q = dict(where="cwa IN (" + ",".join(f"'{o}'" for o in OFFICES) + ")", outFields="cwa,city,state", outSR=4326,
             maxAllowableOffset=0.001, geometryPrecision=5, f="geojson")
    gj = http_get(CWA_URL + "?" + urlencode(q), timeout=(15, 180)).json()
    got = {f["properties"]["cwa"].strip() for f in gj.get("features", [])}
    if got != set(OFFICES): raise RuntimeError(f"NWS reference service returned CWAs {sorted(got)}, expected {sorted(OFFICES)}")
    out.parent.mkdir(parents=True, exist_ok=True); tmp = f"{out}.{os.getpid()}.tmp"
    with open(tmp, "w") as fh: json.dump(gj, fh)
    os.replace(tmp, out)

def ensure_cwa_file():
    f = cwa_file()
    if f is None:
        f = config.WORK / "cwa_new_england.geojson"; print(f"downloading CWA polygons from the NWS reference map service -> {f}", flush=True); fetch_cwa_file(f)
    os.environ["NWSNOW_CWA_FILE"] = str(f)           # worker processes use the same file
    return f

def load_cwa():
    """{office: shapely polygon} for every office in the CWA file"""
    if not _CWA:
        for ft in json.load(open(ensure_cwa_file()))["features"]:
            g = shape(ft["geometry"])
            if not g.is_valid:
                g = make_valid(g)
                if g.geom_type == "GeometryCollection": g = shapely.union_all([p for p in g.geoms if p.geom_type in ("Polygon", "MultiPolygon")])
            _CWA[ft["properties"]["cwa"].strip()] = g
    return _CWA


# ======================================================================================================================================
# STAGE 1 - terrain: AWS Terrain Tiles (Terrarium PNG) regridded to a lon/lat elevation grid (m)  ->  terrain/dem.npz (lon, lat, elev)
# ======================================================================================================================================
DEM_Z = 8
DEM_LON0, DEM_LON1, DEM_LAT0, DEM_LAT1 = -79.6, -63.4, 39.6, 48.7
DEM_STEP = 0.004  # degrees

def tile_xy(lon, lat, z):
    n = 2 ** z
    x = (lon + 180) / 360 * n
    y = (1 - math.asinh(math.tan(math.radians(lat))) / math.pi) / 2 * n
    return x, y

def fetch_tile(t):
    z, x, y = t
    path = str(config.TERRAIN / "tiles" / f"{z}_{x}_{y}.png")
    if not os.path.exists(path):
        os.makedirs(os.path.dirname(path), exist_ok=True)
        url = f"https://s3.amazonaws.com/elevation-tiles-prod/terrarium/{z}/{x}/{y}.png"
        tmp = path + f".{os.getpid()}.part"
        open(tmp, "wb").write(http_get(url).content); os.replace(tmp, path)
    a = np.asarray(Image.open(path).convert("RGB")).astype(np.float32)
    return (x, y), a[..., 0] * 256 + a[..., 1] + a[..., 2] / 256 - 32768

def stage_terrain():
    x0, y1 = tile_xy(DEM_LON0, DEM_LAT0, DEM_Z); x1, y0 = tile_xy(DEM_LON1, DEM_LAT1, DEM_Z)
    xs = range(int(x0), int(x1) + 1); ys = range(int(y0), int(y1) + 1)
    tiles = [(DEM_Z, x, y) for x in xs for y in ys]
    print(len(tiles), "tiles")
    with ThreadPoolExecutor(8) as ex:
        res = dict(ex.map(fetch_tile, tiles))
    W = (len(xs)) * 256; H = (len(ys)) * 256
    mos = np.zeros((H, W), np.float32)
    for (x, y), e in res.items():
        mos[(y - ys[0]) * 256:(y - ys[0] + 1) * 256, (x - xs[0]) * 256:(x - xs[0] + 1) * 256] = e
    lon = np.arange(DEM_LON0, DEM_LON1 + 1e-9, DEM_STEP); lat = np.arange(DEM_LAT1, DEM_LAT0 - 1e-9, -DEM_STEP)
    px = np.array([tile_xy(l, 45, DEM_Z)[0] for l in lon]); px = (px - xs[0]) * 256 - 0.5
    py = np.array([tile_xy(-70, la, DEM_Z)[1] for la in lat]); py = (py - ys[0]) * 256 - 0.5
    # bilinear sample
    ix = np.clip(np.floor(px).astype(int), 0, W - 2); fx = (px - ix)[None, :]
    iy = np.clip(np.floor(py).astype(int), 0, H - 2); fy = (py - iy)[:, None]
    a = mos[iy][:, ix]; b = mos[iy][:, ix + 1]; c = mos[iy + 1][:, ix]; d = mos[iy + 1][:, ix + 1]
    elev = (a * (1 - fx) + b * fx) * (1 - fy) + (c * (1 - fx) + d * fx) * fy
    config.TERRAIN.mkdir(parents=True, exist_ok=True)
    np.savez_compressed(config.DEM, lon=lon, lat=lat, elev=elev.astype(np.float32))
    http_summary()
    print("grid", elev.shape, "min/max", float(elev.min()), float(elev.max()))

# ======================================================================================================================================
# STAGE 2 - cube: NOHRSC sfav2 6-h snowfall (inches) for the New England window
# ======================================================================================================================================
# Source: https://www.nohrsc.noaa.gov/snowfall_v2/data/YYYYMM/sfav2_CONUS_6h_YYYYMMDDHH.tif   (named by END time of the period); one .npy (213 x 350, inches) per 6-h period.
# Fewest-bytes strategy (the file layout changed over the record, so it is probed once per month, not once per file):
#   * LZW-compressed files (about 2022 on, ~0.4 MB each)  -> one whole-file request
#   * uncompressed files (earlier, 5.1 MB each, 1 row per strip) -> one or two multi-range requests for just the New England rows AND
#     columns (~0.3 MB instead of 5.1 MB)
# Every result is range/sanity checked and falls back to a whole-file request if anything looks wrong.  Files known to be absent are
# remembered in nohrsc_cube/_missing.txt.
def _tifffile():
    import tifffile; return tifffile

R0, R1 = config.CUBE_ROWS; C0, C1 = config.CUBE_COLS
ROW_BYTES = 1500 * 4; NR = R1 - R0; NC = C1 - C0; PIX = 4
RANGES_PER_REQUEST = 110          # well under Apache's default limit of 200 ranges per request
_layout = {}; _ml = {}; _gl = threading.Lock(); _missing = set(); _mlock = threading.Lock()

def url_of(s): return f"https://www.nohrsc.noaa.gov/snowfall_v2/data/{s[:6]}/sfav2_CONUS_6h_{s}.tif"

def probe(url):
    """('raw', row_offsets) if rows are uncompressed one-per-strip, else ('whole', None); None if the file does not exist"""
    r = http_get(url, headers={"Range": "bytes=0-9999"}, ok404=True)
    if r is None or len(r.content) < 500: return None
    with _tifffile().TiffFile(io.BytesIO(r.content)) as tf:
        p = tf.pages[0]; comp = p.compression.value if hasattr(p.compression, "value") else int(p.compression)
        offs = p.dataoffsets
        if comp == 1 and p.rowsperstrip == 1 and len(offs) == 850 and offs[1] - offs[0] == ROW_BYTES:
            return ("raw", [int(o) for o in offs[R0:R1]])
    return ("whole", None)

def month_layout(s):
    ym = s[:6]
    with _gl: lock = _ml.setdefault(ym, threading.Lock())
    with lock:
        if ym not in _layout:
            lay = probe(url_of(s))
            if lay is None: return None            # probe file absent: try again with the next file
            _layout[ym] = lay
        return _layout[ym]

def parse_multipart(r, rows_expected):
    """bytes of the requested ranges, concatenated, from a 206 multipart/byteranges response (or a single-range 206)"""
    ct = r.headers.get("Content-Type", "")
    m = re.search(r'boundary=("?)([^";]+)\1', ct)
    if not m:                                               # single range
        return [r.content]
    b = b"--" + m.group(2).encode(); parts = []
    for chunk in r.content.split(b):
        if b"\r\n\r\n" not in chunk: continue
        body = chunk.split(b"\r\n\r\n", 1)[1]
        if body.endswith(b"\r\n"): body = body[:-2]
        if body.endswith(b"--"): body = body[:-2].rstrip(b"\r\n")
        parts.append(body)
    return parts

def fetch_raw(url, offs):
    """New England rows+columns of an uncompressed file with multi-range requests -> float32 (NR, NC) or None"""
    out = []
    for i in range(0, NR, RANGES_PER_REQUEST):
        rows = range(i, min(i + RANGES_PER_REQUEST, NR))
        spec = ",".join(f"{offs[k] + C0 * PIX}-{offs[k] + C1 * PIX - 1}" for k in rows)
        r = http_get(url, headers={"Range": "bytes=" + spec})
        if r.status_code != 206: return None            # server ignored the ranges; use the whole-file path
        parts = parse_multipart(r, len(rows))
        if len(parts) != len(rows) or any(len(p) != NC * PIX for p in parts): return None
        out.extend(parts)
    return np.frombuffer(b"".join(out), dtype="<f4").reshape(NR, NC)

def fetch_whole(url):
    r = http_get(url, ok404=True)
    if r is None: return None
    return _tifffile().imread(io.BytesIO(r.content))[R0:R1, C0:C1]

def sane(a):
    """snowfall in inches: 0..a few hundred, or the -99999 no-data flag.  Anything else means the layout guess was wrong."""
    return bool(np.all((a < -9000) | ((a > -0.01) & (a < 400))))

def fetch_period(t):
    s = t.strftime("%Y%m%d%H"); out = config.CUBE / (s + ".npy")
    if out.exists(): return "have"
    if s in _missing: return "missing"
    try:
        lay = month_layout(s)
        if lay is None: _mark_missing(s); return "missing"
        a = None
        if lay[0] == "raw":
            a = fetch_raw(url_of(s), lay[1])
            if a is not None and not sane(a): a = None
        if a is None: a = fetch_whole(url_of(s))
        if a is None: _mark_missing(s); return "missing"
        a = np.array(a, dtype="f4"); a[a < -9000] = np.nan
        tmp = str(out) + f".{os.getpid()}.{threading.get_ident()}.tmp.npy"; np.save(tmp, a); os.replace(tmp, out); return "ok"
    except Exception as e:
        return "err " + str(e)[:80]

def _mark_missing(s):
    with _mlock:
        _missing.add(s)
        with open(config.CUBE / "_missing.txt", "a") as f: f.write(s + "\n")

def cube_periods():
    t = config.START.replace(hour=6); out = []
    while t <= config.END.replace(hour=18):
        if t.month in config.WINTER_MONTHS: out.append(t)
        t += dt.timedelta(hours=6)
    return out

def stage_cube(workers=6):
    config.CUBE.mkdir(parents=True, exist_ok=True)
    mf = config.CUBE / "_missing.txt"
    if mf.exists(): _missing.update(x.strip() for x in open(mf) if x.strip())
    ts = cube_periods(); todo = [t for t in ts if not (config.CUBE / (t.strftime("%Y%m%d%H") + ".npy")).exists()]
    print(f"{len(ts)} six-hour periods, {ts[0]:%Y-%m-%d} .. {ts[-1]:%Y-%m-%d}; {len(ts) - len(todo)} already cached, {len(todo)} to fetch", flush=True)
    c = Counter(); t0 = time.time()
    with ThreadPoolExecutor(workers) as ex:
        for i, r in enumerate(ex.map(fetch_period, todo)):
            c[r.split()[0]] += 1
            if i % 300 == 299:
                el = time.time() - t0; print(f"  {i+1}/{len(todo)} {dict(c)}  ({el/60:.1f} min, ~{el/(i+1)*(len(todo)-i-1)/60:.0f} min left)", flush=True)
    print("final", dict(c), flush=True); http_summary()
    bad = c["err"]
    if bad: print(f"WARNING: {bad} periods failed; re-run this stage to retry them", flush=True)

# ======================================================================================================================================
# STAGE 3 - storms: detect storms from the NOHRSC cube, fetch the matching WFO-edited NDFD snowfall forecasts, one pickle per office
# ======================================================================================================================================
# Forecast source: NCEI NDFD archive (THREDDS), WMO header YSUZ98 "Snowfall Accumulation", 6-h messages on the CONUS 2.5 km grid
#   https://www.ncei.noaa.gov/thredds/catalog/model-ndfd-file_kwbn[-old]/...   (June 2018 onward)
# Observed: the NOHRSC cube made by stage 2.   Output: storms/ndfd_<OFFICE>.pkl with keys R (per-storm dict), CLON, CLAT, CMASK, sl.
# Efficiency: a 3.8 MB day catalog is fetched only for the day of each target time (the previous day only if no earlier run exists that
# day); only the 12-h-lead run is downloaded up front and the 6-h fallback only for storms it cannot cover; each run file is downloaded
# once and decoded once, for the union of all storms/offices that use it, in parallel processes; the 2.5 km -> 0.04 deg area-average
# regrid is one sparse-matrix multiply per message (weights built once per process).
NDFD_BASE = "https://www.ncei.noaa.gov/thredds/"
CUBE_LON = config.NOHRSC_X0 + config.NOHRSC_STEP * (C0 + np.arange(C1 - C0) + 0.5)      # cell-centre longitudes of the cube
CUBE_LAT = config.NOHRSC_Y0 - config.NOHRSC_STEP * (R0 + np.arange(R1 - R0) + 0.5)      # cell-centre latitudes (north to south)
LON2, LAT2 = np.meshgrid(CUBE_LON, CUBE_LAT)
SIX = dt.timedelta(hours=6)

# ---------------------------------------------------------------- NDFD catalog / file access
def day_catalog(day):
    """[(run datetime UTC, urlPath)] of CONUS 6-h snowfall files (WMO YSUZ98) issued on UTC date `day` (cached as small json)."""
    f = config.CATALOG / f"{day:%Y%m%d}.json"
    if f.exists(): return [(dt.datetime.strptime(a, "%Y%m%d%H%M"), b) for a, b in json.load(open(f))]
    old = day < dt.date(2020, 6, 1)      # archive moved collections in June 2020
    url = NDFD_BASE + ("catalog/model-ndfd-file_kwbn-old/" if old else "catalog/model-ndfd-file_kwbn/access/") + f"{day:%Y%m}/{day:%Y%m%d}/catalog.xml"
    r = http_get(url, ok404=True)
    s = r.text if r is not None else ""
    it = re.findall(r'<dataset name="Snowfall Accumulation run on (\S+) (\d\d):(\d\d) WMO: YSUZ98 CC: KWBN"[^>]*urlPath="([^"]*)"', s)
    out = sorted((dt.datetime.strptime(d + h + m, "%Y-%m-%d%H%M").strftime("%Y%m%d%H%M"), u) for d, h, m, u in it)
    json.dump(out, open(f, "w"))
    return [(dt.datetime.strptime(a, "%Y%m%d%H%M"), b) for a, b in out]

def pick_run(target):
    """latest NDFD snowfall run issued at or before `target`; the previous day's catalog is read only if needed"""
    for d in (target.date(), target.date() - dt.timedelta(days=1)):
        runs = [r for r in day_catalog(d) if r[0] <= target]
        if runs: return max(runs)
    return None

def run_path(u): return config.RUNS / (os.path.basename(u) + ".grib2")

def fetch_run(u):
    f = run_path(u)
    if not f.exists(): http_download(NDFD_BASE + "fileServer/" + u, f)
    return f

# ---------------------------------------------------------------- NDFD -> NOHRSC grid (runs in worker processes)
OFFY = np.array([-0.7, -0.35, 0.0, 0.35, 0.7]); OFFX = np.array([-0.5, -0.25, 0.0, 0.25, 0.5])   # NDFD index offsets spanning one 0.04-deg cell
_W = {}
def _weights(m):
    """sparse matrix W: flattened cropped NDFD grid -> flattened NOHRSC cube grid; each cell is the mean of 5x5 bilinear samples"""
    key = (m.Ni, m.Nj, m.projString)
    if key in _W: return _W[key]
    lat, lon = m.latlons(); t = Transformer.from_crs("EPSG:4326", m.projString, always_xy=True)
    x, y = t.transform(lon[[0, 0], [0, 1]], lat[[0, 0], [0, 1]]); x1, y1 = t.transform(lon[1, 0], lat[1, 0])
    x0, y0 = x[0], y[0]; dx = x[1] - x[0]; dy = y1 - y0
    xs, ys = t.transform(LON2, LAT2); fj = (ys - y0) / dy; fi = (xs - x0) / dx
    j0c = max(int(np.floor(fj.min() - 2)), 0); j1c = min(int(np.ceil(fj.max() + 2)) + 1, m.Nj)
    i0c = max(int(np.floor(fi.min() - 2)), 0); i1c = min(int(np.ceil(fi.max() + 2)) + 1, m.Ni); nci = i1c - i0c
    n = fj.size; rows = np.arange(n); R_, C_, V_ = [], [], []
    for a in OFFY:
        for b in OFFX:
            j = np.clip(fj + a, 0, m.Nj - 1).ravel(); i = np.clip(fi + b, 0, m.Ni - 1).ravel()
            jb = np.minimum(np.floor(j).astype(int), m.Nj - 2); ib = np.minimum(np.floor(i).astype(int), m.Ni - 2); wj = j - jb; wi = i - ib
            for dj, di, w in ((0, 0, (1 - wj) * (1 - wi)), (0, 1, (1 - wj) * wi), (1, 0, wj * (1 - wi)), (1, 1, wj * wi)):
                R_.append(rows); C_.append((jb + dj - j0c) * nci + (ib + di - i0c)); V_.append(w / 25.0)
    W = sp.csr_matrix((np.concatenate(V_).astype("f4"), (np.concatenate(R_), np.concatenate(C_))), shape=(n, (j1c - j0c) * nci))
    _W[key] = (W, (slice(j0c, j1c), slice(i0c, i1c))); return _W[key]

def regrid_run(path, ends):
    """{end datetime: float32 array on the cube grid, inches} for the requested 6-h period end times found in the file"""
    import pygrib
    out = {}; want = set(ends)
    g = pygrib.open(str(path))
    for m in g:
        end = dt.datetime(m.validityDate // 10000, m.validityDate // 100 % 100, m.validityDate % 100, m.validityTime // 100)
        if end in want and end not in out:        # some archive files hold duplicate copies of every message
            W, (sj, si) = _weights(m)
            v = np.ma.filled(m.values, 0.0).astype("f4")[sj, si] * np.float32(39.3701)          # m -> in
            out[end] = (W @ v.ravel()).reshape(LON2.shape).astype("f4")
    g.close(); return out

# ---------------------------------------------------------------- cube and storm detection

def load_cube():
    files = sorted(f for f in glob.glob(str(config.CUBE / "??????????.npy")) if len(os.path.basename(f)) == 14)
    if not files: sys.exit("No NOHRSC files found - run the 'cube' stage first.")
    ts = [dt.datetime.strptime(os.path.basename(f)[:10], "%Y%m%d%H") for f in files]
    cube = np.empty((len(files), R1 - R0, C1 - C0), "f4")          # preallocated: no second copy of the whole cube
    for i, f in enumerate(files): cube[i] = np.load(f)
    return cube, ts

def detect(off, cwa, cube, ts):
    M = contains_xy(cwa[off], LON2, LAT2); sub = cube[:, M]; fin = np.isfinite(sub)
    cnt = fin.sum(1); ser = np.where(cnt > 0.5 * M.sum(), np.where(fin, sub, 0).sum(1) / np.maximum(cnt, 1), np.nan)
    hot = np.where(np.nan_to_num(ser) >= config.THR_SERIES)[0]
    grp = []
    for i in hot:
        if grp and ts[i] - ts[grp[-1][-1]] <= dt.timedelta(hours=config.GAP_HOURS): grp[-1].append(i)
        else: grp.append([i])
    ev = []
    for g in grp:
        i0, i1 = g[0], g[-1]; tot = np.nansum(ser[i0:i1 + 1]); w0 = ts[i0] - SIX; w1 = ts[i1]
        if tot < config.THR_TOTAL: continue
        if (w1 - w0).total_seconds() / 3600 > config.MAX_HOURS:         # split very long events at their quietest 6-h period
            seg = [(i0, i1)]
            while True:
                j0, j1 = seg[-1]
                if (ts[j1] - ts[j0]).total_seconds() / 3600 + 6 <= config.MAX_HOURS: break
                k = j0 + 4 + int(np.nanargmin(ser[j0 + 4:j0 + 8])); seg[-1] = (j0, k); seg.append((k + 1, j1))
            for a, b in seg: ev.append((ts[a] - SIX, ts[b], np.nansum(ser[a:b + 1])))
        else: ev.append((w0, w1, tot))
    print(f"{off}: {len(ev)} storms (CWA-mean total >= {config.THR_TOTAL} in)", flush=True)
    return M, ev

def ends_of(w0, w1): return [w0 + SIX * (i + 1) for i in range(int((w1 - w0).total_seconds() // 21600))]

def process_runs(runs, needs, jobs):
    """download + decode every run file in `runs` for the end times in needs[run]; returns {run url: {end: array}}"""
    runs = sorted(runs)
    with ThreadPoolExecutor(min(len(runs), 8) or 1) as ex: list(ex.map(fetch_run, runs))      # pacing/connection caps live in http_get()
    paths = {u: run_path(u) for u in runs if run_path(u).exists()}
    if jobs > 1 and len(paths) > 1:
        with ProcessPoolExecutor(min(jobs, len(paths))) as ex:
            res = list(ex.map(regrid_run, paths.values(), [sorted(needs[u]) for u in paths], chunksize=1))
    else:
        res = [regrid_run(p, sorted(needs[u])) for u, p in paths.items()]
    return dict(zip(paths, res))

def stage_storms(offices, workers=8, jobs=None):
    jobs = jobs or os.cpu_count() or 1
    for d in (config.CATALOG, config.RUNS, config.STORMS): d.mkdir(parents=True, exist_ok=True)
    cube, ts = load_cube(); idx = {t: i for i, t in enumerate(ts)}
    print("cube", cube.shape, ts[0], "..", ts[-1], flush=True)
    cwa = load_cwa(); plan = {off: detect(off, cwa, cube, ts) for off in offices}
    jobs_list = [(off, w0, w1) for off in offices for w0, w1, tot in plan[off][1] if all(p in idx for p in ends_of(w0, w1))]
    lead_runs = {}
    def assign(lead, items):
        days = sorted({(w0 - dt.timedelta(hours=lead)).date() for _, w0, _ in items})
        todo = [d for d in days if not (config.CATALOG / f"{d:%Y%m%d}.json").exists()]
        print(f"lead {lead} h: {len(todo)} catalog days to fetch ({len(days) - len(todo)} cached)", flush=True)
        with ThreadPoolExecutor(workers) as ex: list(ex.map(day_catalog, todo))
        return {(off, w0, w1): pick_run(w0 - dt.timedelta(hours=lead)) for off, w0, w1 in items}
    totals = {}                                                    # (off, w0, w1) -> (total array, run, lead)
    pending = list(jobs_list)
    for lead in config.LEADS:
        if not pending: break
        pick = {k: v for k, v in assign(lead, pending).items() if v is not None}
        needs = {}
        for (off, w0, w1), (rt, u) in pick.items(): needs.setdefault(u, set()).update(ends_of(w0, w1))
        print(f"lead {lead} h: {len(needs)} NDFD run files for {len(pick)} storms", flush=True)
        dec = process_runs(needs.keys(), needs, jobs)
        nxt = []
        for k in pending:
            if k not in pick: nxt.append(k); continue
            rt, u = pick[k]; off, w0, w1 = k; ends = ends_of(w0, w1); d = dec.get(u, {})
            if all(e in d for e in ends): totals[k] = (sum(d[e] for e in ends), rt, lead)
            else: nxt.append(k)
        pending = nxt
    for off in offices:
        M, ev = plan[off]; res = {}
        for w0, w1, tot in ev:
            k = (off, w0, w1)
            if k not in totals: continue
            fc, rt, lead = totals[k]; truth = np.sum(cube[[idx[p] for p in ends_of(w0, w1)]], axis=0)
            res[w0.strftime("%Y%m%d%H")] = dict(fc=fc, nsa=truth, label=w0.strftime("%d %b %Y").lstrip("0"),
                          issuance="NDFD run %s Z (lead %d h)" % (rt.strftime("%d %b %H:%M"), lead),
                          window="%s – %s" % (w0.strftime("%HZ %d %b %Y"), w1.strftime("%HZ %d %b %Y")),
                          w0=w0, w1=w1, lead=lead, run=rt, cwa_mean_nsa=tot)
        ys, xs = np.where(M); sl = (slice(max(ys.min() - 2, 0), ys.max() + 3), slice(max(xs.min() - 2, 0), xs.max() + 3))
        out = dict(R={k: dict(v, fc=v["fc"][sl], nsa=v["nsa"][sl]) for k, v in res.items()}, CLON=CUBE_LON[sl[1]], CLAT=CUBE_LAT[sl[0]], CMASK=M[sl], sl=sl)
        pickle.dump(out, open(config.STORMS / f"ndfd_{off}.pkl", "wb"))
        print(f"{off}: saved {len(res)} storms ({len(ev) - len(res)} skipped: missing NOHRSC period or no covering NDFD run)", flush=True)
    http_summary()

# ======================================================================================================================================
# MAP HELPERS - projection, terrain shading, cached base map, boundary styles, fast image warp, color scales
# ======================================================================================================================================
PC = ccrs.PlateCarree()
PROJ = ccrs.LambertConformal(central_longitude=-71.5, central_latitude=44.0, standard_parallels=(42, 46))

def load_dem():
    z = np.load(config.DEM)
    return z["lon"], z["lat"], z["elev"]

def hillshade_rgb(lon, lat, elev, strength=1.0, vert_exag=6.0):
    """Hypsometric tint multiplied by hillshade. Returns (rgb, extent)."""
    e = np.where(elev < 0, 0, elev)
    dy = abs(lat[1] - lat[0]) * 111_000.0
    dx = abs(lon[1] - lon[0]) * 111_000.0 * np.cos(np.radians(np.nanmean(lat)))
    ls = LightSource(azdeg=315, altdeg=40)
    sh = ls.hillshade(e, vert_exag=vert_exag, dx=dx, dy=dy)
    cmap = LinearSegmentedColormap.from_list("hyps", [
        (0.00, "#e8efd8"), (0.08, "#dfe6c6"), (0.20, "#d7d6ab"), (0.38, "#c9b98f"),
        (0.60, "#b39d7d"), (0.85, "#9b8a78"), (1.00, "#ffffff")])
    tint = cmap(np.clip(e / 1500.0, 0, 1))[..., :3]
    shade = (0.50 + 0.50 * sh)[..., None]
    rgb = np.clip(tint * (1 - strength + strength * shade * 1.15), 0, 1)
    ext = (lon[0], lon[-1], lat[-1], lat[0])
    return rgb, ext

def draw_terrain(ax, rgb, ext, regrid=1600):
    imshow_proj(ax, rgb, ext, regrid_shape=regrid, interpolation="bilinear", zorder=0)

# ---- fast base map: Natural Earth features are clipped to the study domain and projected ONCE (cached on disk), then drawn as
# plain matplotlib collections in map coordinates.  This avoids cartopy re-projecting the huge 10 m ocean/land polygons on every figure.
_BM = {}
BOX = (-80.5, 38.5, -62.5, 49.5)          # lon0, lat0, lon1, lat1: covers every office map with margin
_FEATS = {"ocean": ("physical", "ocean"), "lakes": ("physical", "lakes"), "lakes_na": ("physical", "lakes_north_america"),
          "coast": ("physical", "coastline"), "states": ("cultural", "admin_1_states_provinces_lines"),
          "borders": ("cultural", "admin_0_boundary_lines_land")}

MIN_ISLET = 5e-4          # deg^2 (~5 km^2): smaller islets are generalized away (outline and fill) so the coast is legible at map scale
MIN_COAST_LEN = 0.03      # deg: shorter open coastline fragments are dropped

def _generalize(g):
    """drop tiny holes (islets in water polygons) and tiny rings / short fragments (coast lines)"""
    import shapely
    from shapely.geometry import Polygon, MultiPolygon, LineString, MultiLineString
    parts = list(g.geoms) if hasattr(g, "geoms") else [g]; out = []
    for p in parts:
        if p.is_empty: continue
        if p.geom_type == "Polygon":
            out.append(Polygon(p.exterior, [r for r in p.interiors if Polygon(r).area >= MIN_ISLET]))
        elif p.geom_type == "LineString":
            if p.is_closed and len(p.coords) > 3 and Polygon(p).area < MIN_ISLET: continue
            if not p.is_closed and p.length < MIN_COAST_LEN: continue
            out.append(p)
    if not out: return None
    if out[0].geom_type == "Polygon": return MultiPolygon(out)
    return MultiLineString(out)

def _build_basemap():
    import shapely, pickle
    box = shapely.box(*BOX); out = {}
    def proj(c): return PROJ.transform_points(PC, c[:, 0], c[:, 1])[:, :2]
    for key, (cat, name) in _FEATS.items():
        geoms = []
        for g in cfeature.NaturalEarthFeature(cat, name, "10m").geometries():
            if g.is_empty or not g.intersects(box): continue
            g = shapely.clip_by_rect(g, *BOX)
            if key in ("ocean", "lakes", "lakes_na", "coast"): g = _generalize(g)
            if g is None or g.is_empty: continue
            geoms.append(shapely.transform(shapely.segmentize(g, 0.02), proj))
        out[key] = geoms
    return out

def basemap():
    if not _BM:
        import pickle
        f = config.TERRAIN / "basemap_lcc_v2.pkl"
        if f.exists(): _BM.update(pickle.load(open(f, "rb")))
        else:
            _BM.update(_build_basemap()); f.parent.mkdir(parents=True, exist_ok=True); tmp = str(f) + f".{__import__('os').getpid()}.tmp"; pickle.dump(dict(_BM), open(tmp, "wb")); __import__('os').replace(tmp, f)
    return _BM

def _polys(geoms):
    from matplotlib.path import Path
    paths = []
    for g in geoms:
        for p in (g.geoms if hasattr(g, "geoms") else [g]):
            if p.geom_type != "Polygon" or p.is_empty: continue
            rings = [np.asarray(p.exterior.coords)] + [np.asarray(r.coords) for r in p.interiors]
            v = np.concatenate(rings); c = np.concatenate([[Path.MOVETO] + [Path.LINETO] * (len(r) - 2) + [Path.CLOSEPOLY] for r in rings])
            paths.append(Path(v, c))
    return paths

def _lines(geoms):
    out = []
    for g in geoms:
        for p in (g.geoms if hasattr(g, "geoms") else [g]):
            if p.geom_type == "LineString" and not p.is_empty: out.append(np.asarray(p.coords))
            elif p.geom_type == "Polygon" and not p.is_empty: out.append(np.asarray(p.exterior.coords))
    return out

# Line hierarchy (cartographic convention: heavier/darker = higher-order boundary), all drawn ABOVE the data layer (zorder 4):
#   international boundary  solid, dark, white casing        state / province  dashed, mid-gray        coastline  thin solid, blue
BOUNDARY_STYLE = dict(
    coast=dict(color="#1f5f8b", lw=0.9, ls="-"),
    states=dict(color="#4a4a4a", lw=0.9, ls=(0, (5, 2.5))),
    borders=dict(color="#111111", lw=1.5, ls="-"))

def draw_base(ax, water_color="#cfe3f1", lw_state=None):
    """water fill under the data (zorder 1); coastline, state/province and international boundaries over the data (zorder 5.1-5.3)"""
    from matplotlib.collections import PathCollection, LineCollection
    bm = basemap()
    water = _polys(bm["ocean"] + bm["lakes"] + bm["lakes_na"])
    ax.add_collection(PathCollection(water, facecolor=water_color, edgecolor="none", zorder=1), autolim=False)
    ax.add_collection(LineCollection(_lines(bm["states"]), colors=BOUNDARY_STYLE["states"]["color"], linewidths=BOUNDARY_STYLE["states"]["lw"],
                                     linestyles=BOUNDARY_STYLE["states"]["ls"], zorder=5.1), autolim=False)
    ax.add_collection(LineCollection(_lines(bm["coast"]), colors=BOUNDARY_STYLE["coast"]["color"], linewidths=BOUNDARY_STYLE["coast"]["lw"], zorder=5.2), autolim=False)
    b = _lines(bm["borders"])
    ax.add_collection(LineCollection(b, colors="white", linewidths=BOUNDARY_STYLE["borders"]["lw"] + 1.8, zorder=5.25), autolim=False)      # casing
    ax.add_collection(LineCollection(b, colors=BOUNDARY_STYLE["borders"]["color"], linewidths=BOUNDARY_STYLE["borders"]["lw"], zorder=5.3), autolim=False)

# ---- fast image warp: lon/lat grid image -> map projection, nearest-neighbour, index map cached per (extent, size)
_WARP = {}
def imshow_proj(ax, img, extent, regrid_shape=1400, **kw):
    """Equivalent of cartopy's imshow(transform=PlateCarree, regrid_shape=N) but ~100x faster after the first call.
    img row 0 is north; extent = (lon0, lon1, lat_bottom, lat_top)."""
    from pyproj import Transformer
    x0, x1 = ax.get_xlim(); y0, y1 = ax.get_ylim(); w, h = x1 - x0, y1 - y0
    nx, ny = (regrid_shape, max(1, int(round(regrid_shape * h / w)))) if w >= h else (max(1, int(round(regrid_shape * w / h))), regrid_shape)
    H, W = img.shape[:2]; key = (tuple(extent), H, W, x0, x1, y0, y1, nx, ny)
    if key not in _WARP:
        X = x0 + (np.arange(nx) + 0.5) * w / nx; Y = y1 - (np.arange(ny) + 0.5) * h / ny
        XX, YY = np.meshgrid(X, Y)
        lon, lat = Transformer.from_crs(PROJ, 4326, always_xy=True).transform(XX, YY)
        ci = np.floor((lon - extent[0]) / (extent[1] - extent[0]) * W).astype(np.int32)
        ri = np.floor((extent[3] - lat) / (extent[3] - extent[2]) * H).astype(np.int32)
        valid = (ci >= 0) & (ci < W) & (ri >= 0) & (ri < H)
        _WARP[key] = (np.clip(ri, 0, H - 1), np.clip(ci, 0, W - 1), valid)
    ri, ci, valid = _WARP[key]
    out = img[ri, ci]
    if out.ndim == 3:
        out = np.asarray(out, dtype=np.float32 if out.dtype.kind == "f" else out.dtype)
        a = np.where(valid, 1.0 if out.dtype.kind == "f" else 255, 0).astype(out.dtype)
        out = np.concatenate([out[..., :3], a[..., None]], axis=2)
    else:
        out = np.where(valid, out, np.nan)
    kw.pop("transform", None)
    im = ax.imshow(out, origin="upper", extent=(x0, x1, y0, y1), transform=PROJ, **kw)
    ax.set_xlim(x0, x1); ax.set_ylim(y0, y1)
    return im

def outline_geom(g, min_area=0.02, tol=0.008):
    from shapely.geometry import MultiPolygon
    parts = list(g.geoms) if hasattr(g, "geoms") else [g]
    keep = [p for p in parts if p.area >= min_area] or parts
    out = MultiPolygon(keep).simplify(tol, preserve_topology=True)
    return out


# ---- colormaps (diverging for bias, sequential for error magnitude) and the shared terrain background
DIV = LinearSegmentedColormap.from_list("div", ["#1b5fa8", "#6fa3d6", "#cfe0f1", "#f0efec", "#f4d3c4", "#e48a6c", "#b8321f"])
SEQ = LinearSegmentedColormap.from_list("seq", ["#f6f1fa", "#d4c1e6", "#a98bd0", "#7a4fb0", "#4a2682", "#2a1252"])
_cache = {}
def terrain_rgb():
    """hillshaded terrain image, computed once per process"""
    if "t" not in _cache:
        dlon, dlat, elev = load_dem()
        _cache["t"] = hillshade_rgb(dlon, dlat, elev, strength=1.0, vert_exag=6.0)
    return _cache["t"]

# ======================================================================================================================================
# CARTOGRAPHY - furniture shared by every map
# ======================================================================================================================================
# Adds to every map panel, in order of the cartographic hierarchy: graticule (labelled inside the neatline), neighbouring NWS CWA boundaries
# (dash-dot, labelled with the office id), state / province / lake / ocean names (placed outside the data area), reference cities, the
# international boundary re-drawn above the CWA outline, a boundary legend, a dual scale bar (mi + km), a north arrow (rotated for the
# Lambert conformal convergence at the map centre) and the data / projection credit.   Usage: cartography(ax, ext, focus="GYX", shown=["GYX"]).
# Everything that depends on the final layout (sizes in inches, which corner is emptiest) is computed when the figure is saved, via a small
# hook on Figure.savefig, so subplots_adjust() / colorbars added later do not matter.  Multi-panel figures get the legend, scale bar, north
# arrow and credit on the first map panel only.
HALO = [pe.withStroke(linewidth=2.6, foreground="white")]
_REF = {}

# ---------------------------------------------------------------- reference layers (cached)
def reference():
    """state/province polygons, lake polygons (named), ocean, and populated places, all clipped to config study box (lon/lat)"""
    if _REF: return _REF
    f = config.TERRAIN / "reference_layers.pkl"
    if f.exists(): _REF.update(pickle.load(open(f, "rb"))); return _REF
    box = shapely.box(*BOX); out = {"states": [], "lakes": [], "ocean": [], "cities": []}
    r = shp.Reader(shp.natural_earth("10m", "cultural", "admin_1_states_provinces"))
    for g, rec in zip(r.geometries(), r.records()):
        if rec.attributes.get("admin") in ("United States of America", "Canada") and g.intersects(box):
            out["states"].append((rec.attributes["name"], shapely.clip_by_rect(g, *BOX)))
    for nm in ("lakes", "lakes_north_america"):
        r = shp.Reader(shp.natural_earth("10m", "physical", nm))
        for g, rec in zip(r.geometries(), r.records()):
            n = rec.attributes.get("name") or rec.attributes.get("NAME")
            if n and g.intersects(box): out["lakes"].append((n, shapely.clip_by_rect(g, *BOX)))
    for g in cfeature.NaturalEarthFeature("physical", "ocean", "10m").geometries():
        if g.intersects(box): out["ocean"].append(shapely.clip_by_rect(g, *BOX))
    r = shp.Reader(shp.natural_earth("10m", "cultural", "populated_places"))
    for rec in r.records():
        a = rec.attributes; lo, la = a["LONGITUDE"], a["LATITUDE"]
        if BOX[0] < lo < BOX[2] and BOX[1] < la < BOX[3] and (a.get("ADM0NAME") in ("United States of America", "Canada")):
            out["cities"].append((a["NAME"], lo, la, float(a.get("POP_MAX") or 0)))
    import os; f.parent.mkdir(parents=True, exist_ok=True); tmp = str(f) + f".{os.getpid()}.tmp"; pickle.dump(out, open(tmp, "wb")); os.replace(tmp, f); _REF.update(out); return _REF

# ---------------------------------------------------------------- helpers
def _proj(g):
    return shapely.transform(g, lambda c: PROJ.transform_points(PC, c[:, 0], c[:, 1])[:, :2])

def _candidates(p, n=11):
    """interior points of polygon p ordered from the most interior (farthest from the edge) outward"""
    c0 = shapely.maximum_inscribed_circle(p, 0.02).coords[0]; b = p.bounds
    gx, gy = np.meshgrid(np.linspace(b[0], b[2], n + 2)[1:-1], np.linspace(b[1], b[3], n + 2)[1:-1]); gx, gy = gx.ravel(), gy.ravel()
    k = shapely.contains_xy(p, gx, gy); gx, gy = gx[k], gy[k]
    d = shapely.distance(shapely.points(gx, gy), p.boundary) if len(gx) else np.array([])
    o = np.argsort(-d); return [c0] + [(gx[i], gy[i]) for i in o]

_CB = {}
def _coast_buffer(ref):
    if "b" not in _CB: _CB["b"] = shapely.union_all([g.boundary for g in ref["ocean"]]).buffer(0.035)
    return _CB["b"]

def _nice_le(v):
    e = 10 ** np.floor(np.log10(v)); c = [m * e for m in (1, 2, 2.5, 5, 10) if m * e <= v * 1.0001]
    return max(c) if c else v

def _lonlat_txt(v, axis):
    return f"{abs(v):g}°{('E' if v >= 0 else 'W') if axis == 'lon' else ('N' if v >= 0 else 'S')}"

class _Placer:
    """greedy label collision avoidance in inches; circles for labels, rectangles for the corner furniture"""
    def __init__(self): self.items = []; self.rects = []
    def ok(self, x, y, r):
        for x0, y0, x1, y1 in self.rects:
            if max(x0 - x, 0, x - x1) ** 2 + max(y0 - y, 0, y - y1) ** 2 < r ** 2: return False
        return all((x - a) ** 2 + (y - b) ** 2 >= (r + q) ** 2 for a, b, q in self.items)
    def add(self, x, y, r): self.items.append((x, y, r))

# ---------------------------------------------------------------- the furniture
def _furniture(ax, ext, focus, shown, extra, fig_first, full, reserve=(), states_over_data=False):
    fig = ax.figure; ax.apply_aspect(); pos = ax.get_position(); Wi = pos.width * fig.get_figwidth(); Hi = pos.height * fig.get_figheight()
    x0, x1 = ax.get_xlim(); y0, y1 = ax.get_ylim(); w = x1 - x0; h = y1 - y0
    fs = 9.0 if full else 7.2                                           # base label size, points
    def inch(x, y): return ((x - x0) / w * Wi, (y - y0) / h * Hi)       # projected -> inches inside the axes
    def inside(x, y, m=0.04): return x0 + m * w < x < x1 - m * w and y0 + m * h < y < y1 - m * h
    cwa = load_cwa(); ref = reference()
    lon0, lon1, lat0, lat1 = ext; extbox = shapely.box(lon0, lat0, lon1, lat1)
    data = shapely.union_all([cwa[o] for o in shown]).buffer(0.03) if shown else shapely.Polygon()
    pl = _Placer(); halo = [pe.withStroke(linewidth=2.8 if full else 2.2, foreground="white")]
    neigh = [o for o in cwa if o not in shown and o != focus and cwa[o].intersects(extbox)]
    if fig_first:
        # ---- legend, scale bar, north arrow: choose the emptiest corners
        col = OFFICES[focus]["color"] if focus else "#444444"
        vis = shapely.box(x0, y0, x1, y1); bm = basemap()
        def present(key): return any(g.intersects(vis) for g in bm[key])          # only list symbols that appear on this map
        H = [Line2D([], [], color=col, lw=2.2, label="NWS CWA boundary" if not focus else f"{focus} CWA boundary")]
        if neigh: H.append(Line2D([], [], color="#3d3d3d", lw=1.0, ls=(0, (6, 2, 1.2, 2)), label="Neighboring CWA boundary"))
        if present("borders"): H.append(Line2D([], [], color=BOUNDARY_STYLE["borders"]["color"], lw=BOUNDARY_STYLE["borders"]["lw"], label="International boundary (U.S.–Canada)"))
        if present("states"): H.append(Line2D([], [], color=BOUNDARY_STYLE["states"]["color"], lw=BOUNDARY_STYLE["states"]["lw"], ls=BOUNDARY_STYLE["states"]["ls"], label="State / province boundary"))
        if present("coast"): H.append(Line2D([], [], color=BOUNDARY_STYLE["coast"]["color"], lw=1.2, label="Coastline"))
        H += list(extra)
        lfs = fs - 1.2; ncol = 2 if Wi / Hi > 1.25 else 1; rows = -(-len(H) // ncol)
        leg_w = (0.075 * 38 * lfs / 8 if full else 0.065 * 36 * lfs / 8) * ncol + 0.5; leg_h = rows * (lfs + 4.2) / 72 + 0.22 + (lfs + 4) / 72
        credit = ["Data: NWS NDFD (NCEI), NOHRSC  ·  Terrain: AWS Terrain Tiles", "Boundaries: Natural Earth  ·  Lambert conformal conic"] if full else []
        L_km = _nice_le(0.22 * w / 1000.0); L_mi = _nice_le(0.22 * w / 1609.344)
        bar_in = max(L_km * 1000 / w * Wi, L_mi * 1609.344 / w * Wi); sb_w = max(bar_in + 0.75, (2.9 if full else 0)); sb_h = 0.61 + 0.13 * len(credit)
        # north-arrow direction (grid convergence at the map centre)
        cx, cy = (x0 + x1) / 2, (y0 + y1) / 2; clo, cla = PC.transform_point(cx, cy, PROJ); p2 = PROJ.transform_point(clo, cla + 0.05, PC)
        phi = np.arctan2(p2[0] - cx, p2[1] - cy); na_w, na_h = 0.45, 0.85
        items = {"legend": (leg_w, leg_h), "scale": (sb_w, sb_h), "north": (na_w, na_h)}
        corners = {"ll": (0.0, 0.0), "lr": (1.0, 0.0), "ul": (0.0, 1.0), "ur": (1.0, 1.0)}
        pdata = _proj(shapely.union_all([cwa[o] for o in shown])) if shown else shapely.Polygon()
        def rect(name, corner):
            wi, hi = items[name]; m = 0.08; cxn, cyn = corners[corner]; my = 0.27 if cyn == 0 else m     # bottom clearance keeps the graticule labels visible
            xa = m if cxn == 0 else Wi - m - wi; ya = my if cyn == 0 else Hi - m - hi
            return xa, ya, wi, hi
        def cost(name, corner):
            xa, ya, wi, hi = rect(name, corner)
            r = shapely.box(x0 + xa / Wi * w, y0 + ya / Hi * h, x0 + (xa + wi) / Wi * w, y0 + (ya + hi) / Hi * h)
            return shapely.intersection(r, pdata).area / r.area
        import itertools
        best = min(itertools.permutations(corners, 3), key=lambda pm: sum(cost(n, c) + 0.002 * i for i, (n, c) in enumerate(zip(items, pm))) + (0.02 if pm[1] not in ("lr", "ll") else 0))
        where = dict(zip(items, best))
        for nm_ in items:
            _xa, _ya, _wi, _hi = rect(nm_, where[nm_]); pl.rects.append((_xa - 0.08, _ya - 0.08, _xa + _wi + 0.08, _ya + _hi + 0.08))

    for (rlo, rla) in reserve:
        rx, ry = PROJ.transform_point(rlo, rla, PC); pl.add(*inch(rx, ry), 0.62)
    # ---- neighbouring CWA boundaries (the displayed offices are drawn by the caller)
    if neigh:
        lines = []
        cb = _coast_buffer(ref)
        for o in neigh:      # same simplification as the focus outline; the part that runs along the coast is left to the coastline symbol
            b = shapely.difference(outline_geom(cwa[o]).boundary, cb)
            if not b.is_empty: lines += _lines([_proj(b)])
        ax.add_collection(LineCollection(lines, colors="white", linewidths=2.6, alpha=0.7, zorder=5.55), autolim=False)
        ax.add_collection(LineCollection(lines, colors="#3d3d3d", linewidths=1.0, linestyles=(0, (6, 2, 1.2, 2)), zorder=5.6), autolim=False)
    # ---- graticule with labels inside the neatline
    span = min(lon1 - lon0, (lat1 - lat0) * 1.3); step = 0.5 if span < 2.5 else (1.0 if span < 5.5 else 2.0)
    mer = np.arange(np.ceil((lon0 - 1) / step) * step, lon1 + 1, step); par = np.arange(np.ceil((lat0 - 1) / step) * step, lat1 + 1, step)
    glines = []; lab_b = []; lab_l = []
    for lo in mer:
        la = np.linspace(lat0 - 2, lat1 + 2, 500); p = PROJ.transform_points(PC, np.full_like(la, lo), la)[:, :2]; glines.append(p)
        k = (p[:, 0] > x0 + 0.05 * w) & (p[:, 0] < x1 - 0.05 * w)
        if k.any(): i = np.argmin(np.where(k, np.abs(p[:, 1] - (y0 + 0.012 * h)), 1e30)); lab_b.append((p[i, 0], p[i, 1], _lonlat_txt(lo, "lon"), abs(p[i, 1] - (y0 + 0.012 * h)) < 0.02 * h))
    for la in par:
        lo = np.linspace(lon0 - 3, lon1 + 3, 500); p = PROJ.transform_points(PC, lo, np.full_like(lo, la))[:, :2]; glines.append(p)
        k = (p[:, 1] > y0 + 0.05 * h) & (p[:, 1] < y1 - 0.05 * h)
        if k.any(): i = np.argmin(np.where(k, np.abs(p[:, 0] - (x0 + 0.010 * w)), 1e30)); lab_l.append((p[i, 0], p[i, 1], _lonlat_txt(la, "lat"), abs(p[i, 0] - (x0 + 0.010 * w)) < 0.02 * w))
    ax.add_collection(LineCollection(glines, colors="#555555", linewidths=0.5, linestyles=(0, (1, 3)), alpha=0.7, zorder=5.05), autolim=False)
    def covered(x, y):      # label would sit under the legend / scale bar / north arrow
        xi, yi = inch(x, y); return any(a - 0.05 <= xi <= c_ + 0.05 and b - 0.05 <= yi <= d_ + 0.05 for a, b, c_, d_ in pl.rects)
    for x, y, t, good in lab_b:
        if good and not covered(x, y): ax.text(x, y, t, ha="center", va="bottom", fontsize=fs - 0.8, color="#333", path_effects=halo, zorder=12, clip_on=True)
    for x, y, t, good in lab_l:
        if good and not covered(x, y): ax.text(x, y, t, ha="left", va="center", fontsize=fs - 0.8, color="#333", path_effects=halo, zorder=12, clip_on=True)
    # ---- name labels placed on the largest piece of each feature that lies OUTSIDE the data area
    def place(geom_ll, text, over_data=False, **kw):
        g = shapely.intersection(geom_ll, extbox) if over_data else shapely.difference(shapely.intersection(geom_ll, extbox), data)
        if g.is_empty: return
        parts = [p for p in (g.geoms if hasattr(g, "geoms") else [g]) if p.geom_type == "Polygon"]
        if not parts: return
        p = max(parts, key=lambda q: q.area)
        if p.area < kw.pop("min_frac", 0.03) * extbox.area: return
        r = 0.055 * len(text) * kw.get("fontsize", fs) / 9 + 0.1
        for c0 in _candidates(p):
            x, y = PROJ.transform_point(c0[0], c0[1], PC)
            if not inside(x, y, 0.06): continue
            xi, yi = inch(x, y)
            if pl.ok(xi, yi, r):
                pl.add(xi, yi, r); ax.text(x, y, text, ha="center", va="center", path_effects=halo, zorder=11, clip_on=True, **kw); return
    # ---- neighbour office ids
    for o in neigh:
        g = shapely.difference(shapely.intersection(cwa[o], extbox), data)
        if g.is_empty or g.area < 0.008 * extbox.area: continue
        parts = [p for p in (g.geoms if hasattr(g, "geoms") else [g]) if p.geom_type == "Polygon"]
        for c0 in _candidates(max(parts, key=lambda q: q.area)):
            x, y = PROJ.transform_point(c0[0], c0[1], PC)
            if not inside(x, y, 0.05): continue
            xi, yi = inch(x, y)
            if pl.ok(xi, yi, 0.3):
                pl.add(xi, yi, 0.3); ax.text(x, y, o, ha="center", va="center", fontsize=fs, weight="bold", color="#3d3d3d", path_effects=halo, zorder=11); break
    st = {}
    for n, g in ref["states"]: st[n] = shapely.union(st[n], g) if n in st else g
    for n, g in sorted(st.items(), key=lambda kv: -shapely.intersection(kv[1], extbox).area):
        place(g, n.upper().replace("É", "E"), over_data=states_over_data, fontsize=fs + 0.3, color="#555555", weight="bold", alpha=0.9, min_frac=0.04)
    oc = shapely.union_all(ref["ocean"]) if ref["ocean"] else None
    if oc is not None: place(oc, "ATLANTIC OCEAN", fontsize=fs + 0.3, color="#3b6f93", style="italic", min_frac=0.06)
    lk = {}
    for n, g in ref["lakes"]: lk[n] = shapely.union(lk[n], g) if n in lk else g
    for n, g in sorted(lk.items(), key=lambda kv: -shapely.intersection(kv[1], extbox).area)[:6]:
        place(g, n, fontsize=fs - 1.0, color="#3b6f93", style="italic", min_frac=0.0015)
    # ---- reference cities (largest first, spaced apart)
    ncity = 8 if full else 4; shown_c = 0
    for nm, lo, la, pop in sorted(ref["cities"], key=lambda c: -c[3]):
        if shown_c >= ncity: break
        x, y = PROJ.transform_point(lo, la, PC)
        if not inside(x, y, 0.05): continue
        xi, yi = inch(x, y); r = 0.065 * len(nm) * (fs - 0.5) / 9 + 0.12
        if not pl.ok(xi + r * 0.6, yi, r) or not pl.ok(xi, yi, 0.12): continue
        pl.add(xi + r * 0.6, yi, r); pl.add(xi, yi, 0.12); shown_c += 1
        ax.plot([x], [y], marker="o", ms=3.8 if full else 3.0, mfc="#111", mec="white", mew=0.8, zorder=11, ls="")
        ax.text(x, y, " " + nm, ha="left", va="center", fontsize=fs - 0.5, color="#111", path_effects=halo, zorder=11, clip_on=True)
    # ---- international boundary line on top of the CWA outline, so a CWA edge that follows the border still shows both
    ax.add_collection(LineCollection(_lines(basemap()["borders"]), colors=BOUNDARY_STYLE["borders"]["color"], linewidths=0.9, zorder=8), autolim=False)
    if not fig_first: return
    # legend
    xl, yl, wl, hl = rect("legend", where["legend"]); cn = where["legend"]
    leg = ax.legend(handles=H, loc={"ll": "lower left", "lr": "lower right", "ul": "upper left", "ur": "upper right"}[cn],
                    bbox_to_anchor=((xl if cn[1] == "l" else xl + wl) / Wi, (yl if cn[0] == "l" else yl + hl) / Hi), bbox_transform=ax.transAxes,
                    fontsize=lfs, framealpha=0.94, edgecolor="#999", borderpad=0.55, handlelength=2.6, labelspacing=0.35, borderaxespad=0.0, ncol=ncol,
                    title="Map layers", title_fontsize=lfs)
    leg.set_zorder(30)
    # scale bar
    xa, ya, wi, hi = rect("scale", where["scale"]); T = ax.transAxes
    ax.add_patch(Rectangle((xa / Wi, ya / Hi), wi / Wi, hi / Hi, transform=T, fc="white", ec="#999", lw=0.8, alpha=0.94, zorder=30))
    m = 0.14; cur = ya + 0.07
    for t in credit[::-1]:
        ax.text((xa + m) / Wi, cur / Hi, t, transform=T, fontsize=6.0, color="#444", va="bottom", ha="left", zorder=32); cur += 0.13
    def bar(L, unit_m, yc, lab_above):
        n = 2 if L >= 2 else 1; Lf = L * unit_m / w                      # bar length as a fraction of the axes width
        xs = [(xa + m) / Wi + Lf * k / n for k in range(n + 1)]
        ax.plot([xs[0], xs[-1]], [yc / Hi] * 2, transform=T, color="black", lw=6.0, solid_capstyle="butt", zorder=31)
        for k in range(n): ax.plot([xs[k], xs[k + 1]], [yc / Hi] * 2, transform=T, color="black" if k % 2 == 0 else "white", lw=3.8, solid_capstyle="butt", zorder=32)
        for k, v in enumerate([0] + [L * (j + 1) / n for j in range(n)]):
            lbl = f"{v:g}" + (f" {'mi' if unit_m > 1500 else 'km'}" if k == n else "")
            ax.text(xs[k], (yc + (0.06 if lab_above else -0.06)) / Hi, lbl, transform=T, fontsize=6.8, ha="center", va="bottom" if lab_above else "top", zorder=32)
    yk = cur + 0.19; ym = yk + 0.10
    bar(L_km, 1000.0, yk, False); bar(L_mi, 1609.344, ym, True)
    # north arrow, rotated to true north at the map centre
    xa, ya, wi, hi = rect("north", where["north"])
    ax.add_patch(Rectangle((xa / Wi, ya / Hi), wi / Wi, hi / Hi, transform=T, fc="white", ec="#999", lw=0.8, alpha=0.9, zorder=30))
    L = 0.46; d = np.array([np.sin(phi), np.cos(phi)]); ctr = np.array([xa + wi / 2, ya + hi * 0.43]); tail = ctr - d * L / 2; head = ctr + d * L / 2; nlab = head + d * 0.15
    ax.annotate("", xy=(head[0] / Wi, head[1] / Hi), xytext=(tail[0] / Wi, tail[1] / Hi), xycoords=T, textcoords=T,
                arrowprops=dict(arrowstyle="simple,head_length=0.9,head_width=0.75,tail_width=0.28", fc="#222", ec="white", lw=0.6), zorder=32)
    ax.text(nlab[0] / Wi, nlab[1] / Hi, "N", transform=T, fontsize=10, weight="bold", ha="center", va="center", zorder=32)

# ---------------------------------------------------------------- public API + savefig hook
def cartography(ax, ext, focus=None, shown=(), extra=(), reserve=(), states_over_data=False):
    """register the furniture for this map axes; it is drawn when the figure is saved.
    reserve = [(lon, lat), ...] label positions the caller draws itself (kept clear of other labels);
    states_over_data = also name states that lie under the data (used for the multi-office regional map)"""
    fig = ax.figure; lst = fig.__dict__.setdefault("_cart", [])
    lst.append((ax, ext, focus, list(shown), extra, reserve, states_over_data))

_orig_savefig = Figure.savefig
def _savefig(self, *a, **k):
    lst = self.__dict__.pop("_cart", [])
    if lst:
        ngeo = sum(isinstance(x, GeoAxes) for x in self.axes)
        for n, (ax, ext, focus, shown, extra, reserve, sod) in enumerate(lst):
            _furniture(ax, ext, focus, shown, extra, n == 0, ngeo == 1, reserve, sod)
    return _orig_savefig(self, *a, **k)
Figure.savefig = _savefig

# ======================================================================================================================================
# PER-OFFICE FIGURES
# ======================================================================================================================================
def draw_office_figures(OFF, STAGES, OUT, SEASON=""):
    """Draw the requested figure stages for one office into directory OUT.  Returns False if the office has no storms in SEASON.
    Per-office loader + native-resolution (unsmoothed) rendering of the 0.04-degree NOHRSC-grid fields.
    STAGES: e (single large mean-error map), a (average error maps), b (close enough / which way it misses), c (charts), d (terrain),
            g (largest storms: forecast | observed | error), m (mean forecast | mean observed | mean error).
    SEASON, e.g. "2025-26": keep only storms starting Jul 1 2025 .. Jun 30 2026."""
    EXTRA=[]          # extra legend handles added below (e.g. the radar-exclusion patch)
    cwa = load_cwa(); RGB, EXT = terrain_rgb()
    NAMES={k:f"NWS {v['name']} ({k})" for k,v in config.OFFICES.items()}
    NAME=NAMES[OFF]
    D=pickle.load(open(config.STORMS/f"ndfd_{OFF}.pkl","rb")); R=D["R"]; CL=D["CLON"]; CA=D["CLAT"]; CM=D["CMASK"]
    if SEASON:
        _y=int(SEASON[:4]); R={k:v for k,v in R.items() if dt.datetime(_y,7,1)<=v["w0"]<dt.datetime(_y+1,7,1)}
        if not R: print(f"{OFF}: no storms in season {SEASON}"); return False
        NAME=f"{NAME}, {_y}–{str(_y+1)[2:]} season"
    KEYS=sorted(R)
    ZOOM=8
    dx=abs(CL[1]-CL[0]); dy=abs(CA[1]-CA[0])
    imext=(CL[0]-dx/2, CL[-1]+dx/2, CA[-1]-dy/2, CA[0]+dy/2)
    b=cwa[OFF].bounds; ext=(b[0]-0.25,b[2]+0.25,b[1]-0.2,b[3]+0.2)
    GEOM=outline_geom(cwa[OFF]); col=OFFICES[OFF]["color"]
    FL=imext[0]+(np.arange(len(CL)*ZOOM)+0.5)*dx/ZOOM; FA=imext[3]-(np.arange(len(CA)*ZOOM)+0.5)*dy/ZOOM
    FX,FY=np.meshgrid(FL,FA); FMASK=shapely.contains_xy(GEOM,FX,FY)
    def aspect():
        p=PROJ.transform_points(PC,np.array([ext[0],ext[1],ext[0],ext[1]]),np.array([ext[2],ext[2],ext[3],ext[3]]))
        return (p[:,1].max()-p[:,1].min())/(p[:,0].max()-p[:,0].min())
    ASP=aspect()
    def prep(a):
        a=np.array(a,float); bad=~(CM&np.isfinite(a))
        if bad.all(): return a
        idx=distance_transform_edt(bad,return_distances=False,return_indices=True); return a[tuple(idx)]
    def smooth_native(a,sigma):
        """mask-aware (normalized-convolution) Gaussian on the native grid"""
        a=np.array(a,float); v=CM&np.isfinite(a)
        if sigma<=0: return np.where(v,a,np.nan)
        num=gaussian_filter(np.where(v,a,0.0),sigma,mode="constant"); den=gaussian_filter(v.astype(float),sigma,mode="constant")
        with np.errstate(all="ignore"): return np.where(v&(den>0.05),num/den,np.nan)
    def smooth_up(a,sigma=0.8,lo=None,hi=None):
        f=prep(smooth_native(a,sigma)); f=zoom(f,ZOOM,order=0,mode="nearest",grid_mode=True)   # nearest: native 0.04° cells shown as-is, no interpolation
        if lo is not None: f=np.clip(f,lo,hi)
        return np.where(FMASK,f,np.nan)
    def frame(ax,lw=2.0):
        ax.set_extent(ext,crs=PC); draw_terrain(ax,RGB,EXT,regrid=1400); draw_base(ax)
        ax.add_geometries([GEOM],PC,facecolor="none",edgecolor="white",linewidth=lw+1.0,zorder=6)
        ax.add_geometries([GEOM],PC,facecolor="none",edgecolor=col,linewidth=lw*0.55,zorder=7)
        cartography(ax,ext,focus=OFF,shown=[OFF],extra=EXTRA)
    def show(ax,fine,cmap,vmin=None,vmax=None,norm=None,alpha=0.93):
        return imshow_proj(ax,fine,imext,regrid_shape=2200,cmap=cmap,norm=norm,vmin=None if norm else vmin,vmax=None if norm else vmax,interpolation="nearest",alpha=alpha,zorder=4)
    # terrain predictors
    def terrain():
        z=np.load(config.DEM); lon,lat,elev=z["lon"],z["lat"],z["elev"]
        if lat[1]<lat[0]: lat=lat[::-1]; elev=elev[::-1]
        ix=(CL[None,:]-lon[0])/(lon[1]-lon[0])+0*CA[:,None]; iy=(CA[:,None]-lat[0])/(lat[1]-lat[0])+0*CL[None,:]
        EL=np.maximum(map_coordinates(elev,[iy,ix],order=1,mode="nearest"),0)
        sea=elev<=0.5; la0=np.searchsorted(lat,max(CA.min()-3,lat[0])); la1=np.searchsorted(lat,min(CA.max()+3,lat[-1])); lo0=np.searchsorted(lon,CL.min()-3); lo1=np.searchsorted(lon,min(CL.max()+3,lon[-1]))
        sub=(slice(la0,la1),slice(lo0,lo1)); ds=distance_transform_edt(~sea[sub])*(abs(lat[1]-lat[0])*111.0)
        DK=map_coordinates(ds,[iy-la0,ix-lo0],order=1,mode="nearest")
        return EL,DK

    # ---- radar beam-blockage exclusion (BTV: KCXX Colchester; NOHRSC/MRMS ray artifacts east and southeast of the radar)
    RADAR=config.RADAR_EXCLUSIONS
    XM=np.zeros(CM.shape,bool); XNOTE=""
    if OFF in RADAR:
        _nm,_la0,_lo0,_a0,_a1,_r0=RADAR[OFF]; _LON,_LAT=np.meshgrid(CL,CA)
        _dy=(_LAT-_la0)*111.0; _dx=(_LON-_lo0)*111.0*np.cos(np.radians(_la0)); _rg=np.hypot(_dx,_dy); _az=(np.degrees(np.arctan2(_dx,_dy))+360)%360
        XM=CM&(_az>=_a0)&(_az<=_a1)&(_rg>_r0); CM=CM&~XM
        XNOTE=f"Gray hatched area excluded: radar beam-blockage artifacts in the NOHRSC analysis (radar {_nm}; azimuth {_a0:.0f}°–{_a1:.0f}°, range > {_r0:.0f} km). Excluded cells are not used in any statistic."
    OUT=str(OUT).rstrip("/")+"/"; os.makedirs(OUT,exist_ok=True)
    NS=len(KEYS)
    F=np.stack([R[k]["fc"] for k in KEYS]); T=np.stack([R[k]["nsa"] for k in KEYS]); E=F-T
    ok=np.isfinite(E)&CM[None]
    tol=np.maximum(2.0,0.30*T); hit=(np.abs(E)<=tol); over=(E>tol); under=(E<-tol)
    okS=ok&(np.maximum(F,T)>=4.0); nS=okS.sum(0).astype(float)          # "heavy-snow cases": cell-storms with >=4 in forecast or observed
    n=ok.sum(0).astype(float); n[n<0.9*NS]=np.nan
    yrs=f"{min(R[k]['w0'] for k in KEYS):%b %Y} – {max(R[k]['w1'] for k in KEYS):%b %Y}"
    PW=6.2; OFFN=NAME
    C_H,C_O,C_U="#4f9d69","#d6604d","#4393c3"
    SUB=f"F = official NWS NDFD snowfall forecast issued ~12 h before onset;  O = NOHRSC snowfall analysis;  error = F − O;  N = {NS} storms, {yrs}."
    DEF="Hit: |F − O| ≤ max(2 in, 30% of O).  High bias: F − O above the tolerance (forecast too high).  Low bias: F − O below −tolerance (forecast too low)."
    WSUB,WBOX=210,170
    _su0=smooth_up
    def smooth_up(a,sigma=0.8,lo=None,hi=None):
        f=_su0(a,sigma,lo,hi)      # nmap fills NaN cells from their nearest valid neighbour; re-blank the excluded sector afterwards
        if XM.any(): f=np.where(zoom(XM.astype(float),ZOOM,order=0,mode="nearest",grid_mode=True)>0.5,np.nan,f)
        return f
    _frame0=frame
    XPATCH=Patch(fc="#9a9a9a",ec="#333",hatch="////",alpha=0.8,label="Excluded: radar beam-blockage artifact")
    if XM.any(): EXTRA.append(XPATCH)
    def frame(ax,lw=2.0):
        _frame0(ax,lw)
        if XM.any():
            imshow_proj(ax,np.where(XM,1.0,np.nan),imext,regrid_shape=1400,cmap=ListedColormap(["#9a9a9a"]),vmin=0,vmax=1,alpha=0.55,interpolation="nearest",zorder=6)
            ax.contourf(CL,CA,XM.astype(float),levels=[0.5,1.5],colors="none",hatches=["////"],transform=PC,zorder=7)
    def wrap(t,w=170): return "\n".join(textwrap.fill(p,w) for p in t.split("\n"))
    def box(fig,text,y=0.01,fs=10.5):
        return   # footer descriptions removed (the maps use the space instead)
    def head(fig,title,sub,y=0.985):
        H=fig.get_figheight(); fig.suptitle(title,fontsize=19,weight="bold",y=1-0.12/H,va="top"); fig.text(0.5,1-0.62/H,wrap(sub,WSUB),ha="center",va="top",fontsize=11,color="#444")
    def cbar(fig,ax,im,lab,ticks=None):
        _a=np.asarray(im.get_array(),float); _lo,_hi=im.get_clim(); _lo_c=bool(np.isfinite(_a).any() and np.nanmin(_a)<_lo-1e-9); _hi_c=bool(np.isfinite(_a).any() and np.nanmax(_a)>_hi+1e-9)
        cb=fig.colorbar(im,ax=ax,orientation="horizontal",fraction=0.045,pad=0.025,aspect=32,shrink=0.95,extend={(0,0):"neither",(1,0):"min",(0,1):"max",(1,1):"both"}[(int(_lo_c),int(_hi_c))]); cb.set_label(lab,fontsize=10.5); cb.ax.tick_params(labelsize=9.5)
        if ticks is not None: cb.set_ticks(ticks)
        return cb
    def mapgrid(nr,nc,extra):
        return plt.subplots(nr,nc,figsize=(PW*nc,PW*ASP*nr+extra),subplot_kw=dict(projection=PROJ),squeeze=False)
    def nice(ax):
        for s in ("top","right"): ax.spines[s].set_visible(False)
    def stack(ax,labels,Hh,Oo,Uu,fs=10.5):
        x=np.arange(len(labels)); ax.bar(x,Hh,color=C_H,width=0.7,label="Hit"); ax.bar(x,Oo,bottom=Hh,color=C_O,width=0.7,label="High bias (forecast too high)"); ax.bar(x,Uu,bottom=Hh+Oo,color=C_U,width=0.7,label="Low bias (forecast too low)")
        for i,h in enumerate(Hh):
            if np.isfinite(h): ax.text(i,h/2,f"{h:.0f}%",ha="center",va="center",color="white",fontsize=fs,weight="bold")
        ax.set_xticks(x); ax.set_xticklabels(labels,fontsize=fs-0.5); ax.set_ylim(0,100); ax.set_ylabel("% of cases",fontsize=11); nice(ax)
    tot_pb=100*E[ok].sum()/T[ok].sum(); tot_miss=np.abs(E[ok]).mean()
    hit_all=100*hit[okS].mean(); hi_all=100*over[okS].mean(); lo_all=100*under[okS].mean()
    dirword="too high" if tot_pb>0 else "too low"

    # ---------------------------------------------------------------- A: average error maps
    for STAGE in STAGES:
        if STAGE=="a":
            with np.errstate(all="ignore"):
                mean=np.nanmean(np.where(ok,E,np.nan),0); mae=np.nanmean(np.where(ok,np.abs(E),np.nan),0)
                sE=np.nansum(np.where(ok,E,0),0); sA=np.nansum(np.where(ok,np.abs(E),0),0); sT=np.nansum(np.where(ok,T,0),0)
                pb=100*sE/sT; nm=100*sA/sT
            for a in (mean,mae,pb,nm): a[np.isnan(n)]=np.nan
            thin=sT<1.0*NS; pb[thin]=np.nan; nm[thin]=np.nan
            P=[("Mean Error (Bias):  mean(F − O)",mean,DIV,-2,2,"mean error, in per storm\nred = high bias (forecast too high);  blue = low bias (forecast too low)"),
               ("Mean Absolute Error (MAE):  mean(|F − O|)",mae,SEQ,0,3,"MAE, in per storm"),
               ("Percent Bias:  100·Σ(F − O) / ΣO",pb,DIV,-50,50,"percent bias, %\nred = high bias (forecast too high);  blue = low bias (forecast too low)"),
               ("Normalized MAE:  100·Σ|F − O| / ΣO",nm,SEQ,0,100,"normalized MAE, %")]
            nc=2; nr=2; PW=7.6; TOPI,BOTI=1.3,1.0; fig,axs=mapgrid(nr,nc,TOPI+BOTI+1.0*(nr-1)); Hf=fig.get_figheight()
            for ax,(t,arr,cm,lo,hi,lab) in zip(axs.ravel(),P):
                frame(ax); im=show(ax,smooth_up(arr,0,-1000,1000),cm,lo,hi); ax.set_title(t,fontsize=13,pad=5); cbar(fig,ax,im,lab)
            head(fig,f"{OFFN}: Snowfall Forecast Error Statistics",SUB)
            box(fig,f"Each cell (~4 km, native NOHRSC grid) is computed over all {NS} storms. Σ sums over storms; F = forecast, O = observed. Percent Bias and Normalized MAE are masked where mean O < 1 in per storm.\n"
                    f"{OFF} area-wide: mean error {np.nanmean(E[ok]):+.2f} in, MAE {tot_miss:.2f} in, percent bias {tot_pb:+.1f}%, normalized MAE {100*np.abs(E[ok]).sum()/T[ok].sum():.1f}%, RMSE {np.sqrt(np.mean(E[ok]**2)):.2f} in.",0.012)
            fig.subplots_adjust(left=0.01,right=0.99,top=1-TOPI/Hf,bottom=BOTI/Hf,wspace=0.04,hspace=1.3/(PW*ASP))
            fig.savefig(OUT+f"{OFF}_1_bias_MAE_maps.png",dpi=125,pad_inches=0.1); plt.close(fig)

        # ---------------------------------------------------------------- B: close enough / which way
        if STAGE=="b":
            def rate(x):
                with np.errstate(all="ignore"): return np.where(nS>=8,100*np.where(CM,(x&okS).sum(0),0)/nS,np.nan)
            H,O,U=rate(hit),rate(over),rate(under); lean=O-U
            reg=np.full(H.shape,4.0); reg[(lean>=20)&(H<60)]=2; reg[(lean<=-20)&(H<60)]=3; reg[H>=60]=1; reg[~np.isfinite(H)]=np.nan
            RC=["#4f9d69","#d6604d","#4393c3","#b9b9c4"]
            RL=["Reliable: hit rate ≥ 60%","High bias: hit rate < 60%, high − low ≥ 20 pts","Low bias: hit rate < 60%, low − high ≥ 20 pts","Mixed: hit rate < 60%, |high − low| < 20 pts"]
            PW=7.0; TOPI,BOTI=1.3,1.6; fig,axs=mapgrid(1,3,TOPI+BOTI); Hf=fig.get_figheight()
            ax=axs[0,0]; frame(ax); im=show(ax,smooth_up(H,0,0,100),"Greens",0,100); ax.set_title("Hit Rate (heavy-snow cases)",fontsize=13.5,pad=6)
            cbar(fig,ax,im,"hit rate, %  of cases with ≥4 in forecast or observed\nhit: |F − O| ≤ max(2 in, 30% of O)")
            ax=axs[0,1]; frame(ax); im=show(ax,smooth_up(lean,0,-100,100),DIV,-40,40); ax.set_title("Forecast Tendency:  high-bias rate − low-bias rate",fontsize=13.5,pad=6)
            cbar(fig,ax,im,"percentage points\nred = high bias dominates;  blue = low bias dominates")
            ax=axs[0,2]; frame(ax); imshow_proj(ax,smooth_up(reg,0,0,5),imext,regrid_shape=2200,cmap=ListedColormap(RC),vmin=0.5,vmax=4.5,interpolation="nearest",alpha=0.93,zorder=4)
            ax.set_title("Reliability Category (from hit rate and tendency)",fontsize=13.5,pad=6)
            _d=cbar(fig,ax,im,"."); _d.ax.set_visible(False)   # keeps panel 3 the same size as panels 1-2
            ax.legend(handles=[Patch(fc=c,ec="#444",label=l) for c,l in zip(RC,RL)]+([XPATCH] if XM.any() else []),loc="upper center",bbox_to_anchor=(0.5,-0.03),fontsize=10,frameon=False,ncol=1)
            head(fig,f"{OFFN}: Snowfall Forecast Hit Rate, Tendency and Reliability Category",SUB)
            box(fig,f"{DEF}  Only heavy-snow cases (cell-storms with F ≥ 4 in or O ≥ 4 in) are counted; cells with fewer than 8 such cases are blank.\n"
                    f"{OFF} area-wide, heavy-snow cases: hit rate {hit_all:.0f}%, high-bias rate {hi_all:.0f}%, low-bias rate {lo_all:.0f}%.",0.012)
            fig.subplots_adjust(left=0.01,right=0.99,top=1-TOPI/Hf,bottom=BOTI/Hf,wspace=0.04)
            fig.savefig(OUT+f"{OFF}_2_hit_rate_maps.png",dpi=125,pad_inches=0.1); plt.close(fig)

        # ---------------------------------------------------------------- C: charts
        if STAGE in ("c","d"): EL,DK=terrain()
        if STAGE=="c":
            strata={"Coast\n(<15 km from ocean)":CM&(DK<15),"Inland lowland\n(<200 m)":CM&(DK>=15)&(EL<200),"Hills\n(200–600 m)":CM&(DK>=15)&(EL>=200)&(EL<600),"Mountains\n(≥600 m)":CM&(DK>=15)&(EL>=600)}
            strata={k:v for k,v in strata.items() if v.sum()>=15}
            fig,axs=plt.subplots(2,2,figsize=(17,13.2)); fig.subplots_adjust(left=0.07,right=0.97,top=0.80,bottom=0.20 if XNOTE else 0.17,wspace=0.2,hspace=0.5)
            Tf=T[ok]; Ef=E[ok]; Ff=F[ok]; bins=[0,1,2,4,6,9,12,100]; bl=["<1","1–2","2–4","4–6","6–9","9–12","12+"]
            sel=lambda i:(Tf>=bins[i])&(Tf<bins[i+1])
            # (a) accuracy by amount
            ax=axs[0,0]; hh=[100*hit[ok][sel(i)].mean() for i in range(len(bl))]; oo=[100*over[ok][sel(i)].mean() for i in range(len(bl))]; uu=[100*under[ok][sel(i)].mean() for i in range(len(bl))]
            stack(ax,[f"{b}\n({sel(i).sum()/1000:.0f}k)" for i,b in enumerate(bl)],np.array(hh),np.array(oo),np.array(uu)); ax.set_xlabel("observed snowfall O (in)   (k = thousands of cell-storm cases)",fontsize=10.5)
            ax.set_title("A. Hit / high-bias / low-bias rate\nby observed snowfall (O) bin",fontsize=13,loc="left")
            # (b) terrain
            ax=axs[0,1]; h2,o2,u2=[],[],[]
            for sn,sm in strata.items():
                m=okS&sm[None]; h2.append(100*hit[m].mean()); o2.append(100*over[m].mean()); u2.append(100*under[m].mean())
            stack(ax,list(strata),np.array(h2),np.array(o2),np.array(u2)); ax.set_title("B. Hit / high-bias / low-bias rate\nby terrain class (heavy-snow cases only)",fontsize=13,loc="left")
            # (c) per-storm scatter
            ax=axs[1,0]; fo=np.array([np.nanmean(F[i][ok[i]]) for i in range(NS)]); ob=np.array([np.nanmean(T[i][ok[i]]) for i in range(NS)])
            mx=max(fo.max(),ob.max())*1.05; xx=np.linspace(0,mx,50)
            ax.fill_between(xx,np.minimum(xx*0.75,xx-1),np.maximum(xx*1.25,xx+1),color="#e9efe9",label="|F − O| ≤ max(25% of O, 1 in)")
            ax.plot(xx,xx,color="#222",lw=1.3,label="perfect forecast (F = O)"); ax.scatter(ob,fo,s=34,color="#3b5b8c",alpha=0.75,edgecolor="white",lw=0.6,label="one storm")
            miss=np.argsort(-np.abs(fo-ob))[:4]
            for i in miss: ax.annotate(R[KEYS[i]]["label"],(ob[i],fo[i]),xytext=(6,6 if fo[i]>ob[i] else -12),textcoords="offset points",fontsize=9.5,color="#333")
            ax.set_xlim(0,mx); ax.set_ylim(0,mx); ax.set_aspect("equal"); ax.set_xlabel(f"O: observed snowfall, {OFF}-area mean (in)",fontsize=10.5); ax.set_ylabel(f"F: forecast snowfall, {OFF}-area mean (in)",fontsize=10.5)
            ax.legend(loc="upper left",frameon=False,fontsize=10); nice(ax); inside=100*np.mean(np.abs(fo-ob)<=np.maximum(0.25*ob,1.0))
            ax.set_title(f"C. Storm-mean F vs O ({OFF}-area average, one dot per storm)\n{inside:.0f}% of storms within tolerance band",fontsize=13,loc="left")
            # (d) cell-level mean forecast by observed amount
            ax=axs[1,1]; b=np.arange(0,24,1.0); xm=[];q1=[];md=[];q3=[]
            for lo in b:
                s=(Tf>=lo)&(Tf<lo+1)
                if s.sum()>60: xm.append(lo+0.5); q1.append(np.percentile(Ff[s],25)); md.append(Ff[s].mean()); q3.append(np.percentile(Ff[s],75))
            ax.fill_between(xm,q1,q3,color="#d9d9d9",label="interquartile range of F"); ax.plot(xm,md,color=C_O,lw=2.6,label="mean F"); ax.plot([0,24],[0,24],color="#222",lw=1.3,ls="--",label="perfect forecast (F = O)")
            ax.set_xlim(0,24); ax.set_ylim(0,24); ax.set_aspect("equal"); ax.set_xlabel("O: observed snowfall at grid cell (in)",fontsize=10.5); ax.set_ylabel("F: forecast snowfall at grid cell (in)",fontsize=10.5); ax.legend(loc="upper left",frameon=False,fontsize=10); nice(ax)
            hi_i=[i for i,x in enumerate(xm) if x>=12]; msg=("mean forecast below observed for O ≥ 12 in (low bias in heavy snow)" if hi_i and np.mean([md[i]-xm[i] for i in hi_i])<-1 else "mean forecast tracks observed")
            ax.set_title(f"D. Grid-cell F vs O, mean F by 1-in O bin\n{msg}",fontsize=13,loc="left")
            head(fig,f"{OFFN}: Snowfall Forecast Verification — Hit Rate and F vs O",SUB,0.985)
            fig.legend(handles=[Patch(fc=C_H,label="Hit"),Patch(fc=C_O,label="High bias (forecast too high)"),Patch(fc=C_U,label="Low bias (forecast too low)")],loc="upper center",bbox_to_anchor=(0.5,0.925),ncol=3,frameon=False,fontsize=11.5,title="Categories in A and B",title_fontsize=10.5)
            box(fig,f"{DEF}\nA: all cell-storm cases. B: heavy-snow cases (F ≥ 4 in or O ≥ 4 in). C: each dot averages F and O over the CWA for one storm. D: grid-cell cases binned by O.",0.012)
            fig.savefig(OUT+f"{OFF}_3_error_charts.png",dpi=120,pad_inches=0.1); plt.close(fig)

        # ---------------------------------------------------------------- D: terrain
        if STAGE=="d":
            rng=np.random.default_rng(1); cm=np.where(CM.ravel())[0]; el=EL.ravel()[cm]; dk=DK.ravel()[cm]
            Xk=np.c_[el/130.0,np.log10(1+dk)/0.22]; K=np.exp(-0.5*cdist(Xk,Xk,"sqeuclidean")).astype("f4")
            def sums(idx):
                e=np.where(ok[idx],E[idx],0); o=np.where(ok[idx],T[idx],0); return e.sum(0).ravel()[cm],o.sum(0).ravel()[cm]
            def pbias(idx):
                se,so=sums(idx)
                with np.errstate(all="ignore"): return 100*se/so,so
            def tofit(idx):
                se,so=sums(idx)
                with np.errstate(all="ignore"): return 100*(K@se)/(K@so)
            def tomap(v): a=np.full(CM.shape,np.nan); a.ravel()[cm]=v; return a
            allidx=np.arange(NS); actual,so_all=pbias(allidx); pred=tofit(allidx); good=np.isfinite(actual)&np.isfinite(pred)&(so_all>=NS*1.0)
            r_ter=np.corrcoef(actual[good],pred[good])[0,1]; expl=100*r_ter**2
            rs=[]
            for _ in range(100):
                p=rng.permutation(NS); a,b=pbias(p[:NS//2]),pbias(p[NS//2:]); m=np.isfinite(a[0])&np.isfinite(b[0])&(a[1]>5)&(b[1]>5); rs.append(np.corrcoef(a[0][m],b[0][m])[0,1])
            rmed=np.median(rs); p=np.random.default_rng(7).permutation(NS); a,b=pbias(p[:NS//2]),pbias(p[NS//2:]); mh=np.isfinite(a[0])&np.isfinite(b[0])&(a[1]>5)&(b[1]>5)
            fig=plt.figure(figsize=(21,8.5+6.4*ASP)); gs=fig.add_gridspec(2,3,hspace=0.3,wspace=0.16,height_ratios=[ASP*1.0+0.28,0.9])
            mapspecs=[("A. Percent Bias by grid cell:  100·Σ(F − O) / ΣO",actual),("B. Kernel-regression fit of Percent Bias\non elevation and distance to ocean",pred)]
            for j,(t,v) in enumerate(mapspecs):
                ax=fig.add_subplot(gs[0,j],projection=PROJ); frame(ax); im=show(ax,smooth_up(tomap(np.where(so_all>=NS*1.0,v,np.nan)),0,-1000,1000),DIV,-50,50); ax.set_title(t,fontsize=13,pad=6)
                cbar(fig,ax,im,"percent bias, %\nred = high bias (forecast too high);  blue = low bias (forecast too low)")
            ax=fig.add_subplot(gs[0,2]); ax.hexbin(a[0][mh],b[0][mh],gridsize=45,extent=(-80,80,-80,80),bins="log",cmap="Greys",mincnt=1,linewidths=0); ax.plot([-80,80],[-80,80],color="#222",lw=1.2,ls="--")
            ax.set_xlim(-80,80); ax.set_ylim(-80,80); ax.set_aspect("equal"); ax.axhline(0,color="#999",lw=0.8); ax.axvline(0,color="#999",lw=0.8); nice(ax)
            ax.set_xlabel("Percent Bias, random half A of storms (%)",fontsize=10.5); ax.set_ylabel("Percent Bias, other half B of storms (%)",fontsize=10.5)
            ax.set_title(f"C. Split-half reproducibility of Percent Bias\ngrid-cell scatter, median r = {rmed:.2f} over 100 random splits",fontsize=13,loc="left")
            def boots(var,edges,ax,title,xl):
                xs=[]
                for lo,hi in zip(edges[:-1],edges[1:]):
                    sm=CM&(var>=lo)&(var<hi)
                    if sm.sum()<15: continue
                    per=np.array([[E[i][ok[i]&sm].sum(),np.abs(E[i][ok[i]&sm]).sum(),T[i][ok[i]&sm].sum()] if (ok[i]&sm).sum()>=5 else [np.nan]*3 for i in range(NS)])
                    st=lambda ix:(100*np.nansum(per[ix,0])/np.nansum(per[ix,2]),100*np.nansum(per[ix,1])/np.nansum(per[ix,2]))
                    f=st(allidx); bs=np.array([st(rng.integers(0,NS,NS)) for _ in range(500)]); ci=np.nanpercentile(bs,[5,95],axis=0); xs.append((f,ci,lo,hi))
                xp=np.arange(len(xs))
                for k,(c,lab,off) in enumerate((("#d6604d","Percent Bias",-0.12),("#5e3c99","Normalized MAE",0.12))):
                    y=[x[0][k] for x in xs]; ax.errorbar(xp+off,y,yerr=[[x[0][k]-x[1][0,k] for x in xs],[x[1][1,k]-x[0][k] for x in xs]],fmt="o-",color=c,lw=2.2,ms=7,capsize=4,label=lab)
                ax.axhline(0,color="#444",lw=1); ax.set_xticks(xp); ax.set_xticklabels([("%g–%g"%(x[2],x[3]) if x[3]<1e3 else "%g+"%x[2]) for x in xs],fontsize=11)
                ax.set_xlabel(xl,fontsize=11); ax.set_ylabel("% of observed snowfall (ΣO)",fontsize=11); ax.set_title(title,fontsize=13,loc="left"); ax.legend(frameon=False,fontsize=10.5); nice(ax)
            boots(EL,[-5,50,150,300,500,800,3000],fig.add_subplot(gs[1,0]),"D. Percent Bias and Normalized MAE by elevation band\n(markers = all storms; bars = 90% storm-bootstrap interval)","elevation (m)")
            boots(DK,[0,5,15,30,60,90,1000],fig.add_subplot(gs[1,1]),"E. Percent Bias and Normalized MAE by distance to ocean\n(markers = all storms; bars = 90% storm-bootstrap interval)","distance from the ocean (km)")
            ax=fig.add_subplot(gs[1,2]); ax.axis("off")
            ax.text(0,1,wrap(f"Notes\n\n• Percent Bias = 100·Σ(F − O)/ΣO per grid cell over {NS} storms.\n\n• C: Percent Bias computed separately from two random halves of the storms; high correlation means the spatial pattern is reproducible rather than sampling noise.\n\n• B: kernel regression of Percent Bias on elevation and log distance to ocean; its correlation with A is r = {r_ter:.2f} (r² = {expl:.0f}%), the fraction of the spatial pattern explained by terrain and coast distance alone.\n\n• D, E: a Percent Bias interval that excludes zero indicates bias beyond sampling uncertainty."+(f"\n\n• {XNOTE}" if XNOTE else ""),58),va="top",fontsize=11.5,linespacing=1.4,bbox=dict(boxstyle="round,pad=0.7",fc="#f4f4f1",ec="#cfcfc8"))
            head(fig,f"{OFFN}: Percent Bias vs Elevation and Distance to Ocean",SUB,0.985)
            fig.savefig(OUT+f"{OFF}_4_bias_terrain_coast.png",dpi=115,pad_inches=0.1,bbox_inches="tight"); plt.close(fig)
            print(OFF,"terrain r=%.2f explains %.0f%%; split-half r=%.2f"%(r_ter,expl,rmed))

        # ---------------------------------------------------------------- G: largest storms, three columns
        if STAGE=="g":
            SC=["#bdd6e5","#6caed3","#3281ba","#08519a","#092593","#fdfe96","#fbc305","#fc8706","#db1400","#9e0000","#680001"]
            SL=["0.1–1","1–2","2–3","3–4","4–6","6–8","8–12","12–18","18–24","24–30","30+ in"]
            cmS=ListedColormap(SC); nS_=BoundaryNorm([0.05,1,2,3,4,6,8,12,18,24,30,60],11)
            order=np.argsort([-np.nanmean(T[i][ok[i]]) for i in range(NS)])[:4]
            PW=6.6; TOPI=2.1; XB=0.0; fig,axs=mapgrid(4,3,TOPI+0.2+XB); Hf=fig.get_figheight()
            for r,i in enumerate(order):
                k=KEYS[i]; m=ok[i]; fc=np.where(CM,F[i],np.nan); ob=np.where(CM,T[i],np.nan)
                for c,(arr,t) in enumerate([(fc,f"F: NDFD forecast (run {R[k]['run']:%d %b %H}Z)"),(ob,"O: NOHRSC observed"),(E[i],f"Error F − O:  mean error {E[i][m].mean():+.1f} in, MAE {np.abs(E[i][m]).mean():.1f} in, RMSE {np.sqrt((E[i][m]**2).mean()):.1f} in")]):
                    ax=axs[r,c]; frame(ax,1.8)
                    if c<2: im2=show(ax,np.where(smooth_up(arr,0,0,80)<0.05,np.nan,smooth_up(arr,0,0,80)),cmS,norm=nS_)
                    else: im=show(ax,smooth_up(arr,0,-30,30),DIV,-10,10)
                    ax.set_title(t,fontsize=11.5,pad=4)
                axs[r,0].text(0.02,0.98,f"{R[k]['label']}\n{R[k]['window']}",transform=axs[r,0].transAxes,va="top",fontsize=10.5,weight="bold",bbox=dict(fc="white",ec="#999",alpha=0.9,boxstyle="round,pad=0.3"),zorder=20)
            fig.legend(handles=[Patch(fc=c,ec="#555",label=l) for c,l in zip(SC,SL)],loc="upper left",ncol=11,fontsize=10,frameon=False,bbox_to_anchor=(0.03,1-1.0/Hf),title="Storm-total snowfall, in (columns 1 and 2)",title_fontsize=10.5)
            cax=fig.add_axes([0.70,1-1.75/Hf,0.25,0.006]); cb=fig.colorbar(im,cax=cax,orientation="horizontal"); cb.ax.xaxis.set_ticks_position("top"); cb.ax.xaxis.set_label_position("top"); cb.set_label("column 3: error F − O (in)\nred = high bias, blue = low bias",fontsize=10); cb.ax.tick_params(labelsize=9)
            head(fig,f"{OFFN}: Four Largest Storms — Forecast (F), Observed (O) and Error (F − O)",SUB)
            fig.subplots_adjust(left=0.005,right=0.995,top=1-TOPI/Hf,bottom=(0.2+XB)/Hf,wspace=0.02,hspace=0.14)
            if XNOTE: box(fig,XNOTE,0.012)
            fig.savefig(OUT+f"{OFF}_5_largest_storms_F_O_error.png",dpi=110,pad_inches=0.1); plt.close(fig)

        # ---------------------------------------------------------------- M: mean forecast, mean observed, mean error (actual snowfall amounts)
        if STAGE=="m":
            with np.errstate(all="ignore"):
                mF=np.nanmean(np.where(ok,F,np.nan),0); mO=np.nanmean(np.where(ok,T,np.nan),0); mE=np.nanmean(np.where(ok,E,np.nan),0)
            for a in (mF,mO,mE): a[np.isnan(n)]=np.nan
            vmax=float(np.ceil(np.nanpercentile(np.r_[mF[np.isfinite(mF)],mO[np.isfinite(mO)]],99.5)/2)*2)
            P=[("Mean Forecast Snowfall:  mean(F)",mF,"GnBu",0,vmax,"mean forecast F, in per storm"),
               ("Mean Observed Snowfall:  mean(O)",mO,"GnBu",0,vmax,"mean observed O, in per storm"),
               ("Mean Error (Bias):  mean(F − O)",mE,DIV,-2,2,"mean error, in per storm\nred = high bias (forecast too high);  blue = low bias (forecast too low)")]
            PW=7.0; TOPI,BOTI=1.3,1.0; fig,axs=mapgrid(1,3,TOPI+BOTI); Hf=fig.get_figheight()
            for ax,(t,arr,cm_,lo,hi,lab) in zip(axs.ravel(),P):
                frame(ax); im=show(ax,smooth_up(arr,0,-1000,1000),cm_,lo,hi); ax.set_title(t,fontsize=13.5,pad=6); cbar(fig,ax,im,lab)
            head(fig,f"{OFFN}: Mean Forecast Snowfall, Mean Observed Snowfall and Mean Error",SUB)
            box(fig,f"Each cell (~4 km, native NOHRSC grid) shows the average storm-total snowfall over all {NS} storms: mean(F) forecast, mean(O) observed, and mean(F − O) their difference. Panels 1 and 2 share one color scale.\n"
                    f"{OFF} area-wide averages per storm: mean forecast {np.nanmean(F[ok]):.2f} in, mean observed {np.nanmean(T[ok]):.2f} in, mean error {np.nanmean(E[ok]):+.2f} in.",0.012)
            fig.subplots_adjust(left=0.01,right=0.99,top=1-TOPI/Hf,bottom=BOTI/Hf,wspace=0.04)
            fig.savefig(OUT+f"{OFF}_6_mean_F_O_error_maps.png",dpi=125,pad_inches=0.1); plt.close(fig)

        # ---------------------------------------------------------------- E: one large Mean Error map per CWA
        if STAGE=="e":
            with np.errstate(all="ignore"): mE=np.nanmean(np.where(ok,E,np.nan),0)
            mE[np.isnan(n)]=np.nan; PW=11.0; WSUB,WBOX=130,130
            TOPI,BOTI=1.6,1.0; fig,axs=mapgrid(1,1,TOPI+BOTI); Hf=fig.get_figheight(); ax=axs[0,0]
            frame(ax); im=show(ax,smooth_up(mE,0,-1000,1000),DIV,-2,2); ax.set_title("Mean Error (Bias):  mean(F − O),  inches of snow per storm",fontsize=14,pad=7)
            cbar(fig,ax,im,"mean error, in per storm\nred = high bias (forecast too high);  blue = low bias (forecast too low)")
            v=mE[np.isfinite(mE)&CM]; hi5=100*np.mean(v>0.5); lo5=100*np.mean(v<-0.5); lowmin=f"{np.nanmin(v):+.1f}".replace("-","−")
            head(fig,f"{OFFN}: Mean Snowfall Forecast Error",SUB)
            box(fig,f"Each cell (~4 km, native NOHRSC grid) is the average of F − O over all {NS} storms.  Color scale is the same for every office (−2 to +2 in).\n"
                    f"{OFF} area-wide mean error {np.nanmean(E[ok]):+.2f} in per storm.  {hi5:.0f}% of the CWA has a high bias above +0.5 in; {lo5:.0f}% has a low bias below −0.5 in.  Largest cell-mean errors: {np.nanmax(v):+.1f} in (high) and {lowmin} in (low).",0.012)
            fig.subplots_adjust(left=0.01,right=0.99,top=1-TOPI/Hf,bottom=BOTI/Hf)
            fig.savefig(OUT+f"{OFF}_0_mean_error_map.png",dpi=140,pad_inches=0.1); plt.close(fig)
        print("ok",OFF,STAGE)
        import gc; plt.close('all'); gc.collect()
    return True

# ======================================================================================================================================
# REGIONAL MAP
# ======================================================================================================================================
def office_field(off, season=None, zoom_=8):
    D = pickle.load(open(config.STORMS / f"ndfd_{off}.pkl", "rb")); R = D["R"]; CL, CA, CM = D["CLON"], D["CLAT"], D["CMASK"].copy()
    if season:
        y = int(season[:4]); R = {k: v for k, v in R.items() if dt.datetime(y, 7, 1) <= v["w0"] < dt.datetime(y + 1, 7, 1)}
    if not R: return None
    XM = np.zeros(CM.shape, bool); note = None
    if off in config.RADAR_EXCLUSIONS:
        nm, la0, lo0, a0, a1, r0 = config.RADAR_EXCLUSIONS[off]; LON, LAT = np.meshgrid(CL, CA)
        dy = (LAT - la0) * 111.0; dx = (LON - lo0) * 111.0 * np.cos(np.radians(la0)); az = (np.degrees(np.arctan2(dx, dy)) + 360) % 360
        XM = CM & (az >= a0) & (az <= a1) & (np.hypot(dx, dy) > r0); CM = CM & ~XM
    keys = sorted(R); NS = len(keys)
    E = np.stack([R[k]["fc"] - R[k]["nsa"] for k in keys]); ok = np.isfinite(E) & CM[None]
    n = ok.sum(0).astype(float); n[n < 0.9 * NS] = np.nan
    with np.errstate(all="ignore"): mE = np.nanmean(np.where(ok, E, np.nan), 0)
    mE[np.isnan(n)] = np.nan; area = float(E[ok].mean())
    v = CM & np.isfinite(mE); f = np.where(v, mE, np.nan); bad = ~v
    if not bad.all(): f = f[tuple(distance_transform_edt(bad, return_distances=False, return_indices=True))]
    f = zoom(np.array(f, float), zoom_, order=0, mode="nearest", grid_mode=True)
    dx_ = abs(CL[1] - CL[0]); dy_ = abs(CA[1] - CA[0]); imext = (CL[0] - dx_ / 2, CL[-1] + dx_ / 2, CA[-1] - dy_ / 2, CA[0] + dy_ / 2)
    FL = imext[0] + (np.arange(len(CL) * zoom_) + 0.5) * dx_ / zoom_; FA = imext[3] - (np.arange(len(CA) * zoom_) + 0.5) * dy_ / zoom_
    FX, FY = np.meshgrid(FL, FA); G = outline_geom(load_cwa()[off])
    f = np.where(shapely.contains_xy(G, FX, FY), f, np.nan)
    if XM.any(): f = np.where(zoom(XM.astype(float), zoom_, order=0, mode="nearest", grid_mode=True) > 0.5, np.nan, f)
    return dict(f=f, imext=imext, XM=XM, CL=CL, CA=CA, NS=NS, area=area, geom=G, dates=(min(R[k]["w0"] for k in keys), max(R[k]["w1"] for k in keys)))

def regional_map(offices=None, season_id=None, out_file=None):
    """One regional map: mean snowfall forecast error (F - O) for several offices on a common color scale (cells computed exactly as in figure
    stage e; radar beam-blockage sectors hatched and excluded)."""
    offices = offices or [o for o in config.ORDER if o not in config.REGIONAL_EXCLUDE]
    cwa = load_cwa(); F = {}
    for o in offices:
        r = office_field(o, season_id)
        if r: F[o] = r
        else: print(f"{o}: no storms, skipped")
    offs = list(F); b = shapely.union_all([cwa[o] for o in offs]).bounds; ext = (b[0] - 0.4, b[2] + 0.4, b[1] - 0.35, b[3] + 0.35)
    p = PROJ.transform_points(PC, np.array([ext[0], ext[1], ext[0], ext[1]]), np.array([ext[2], ext[2], ext[3], ext[3]]))
    asp = (p[:, 1].max() - p[:, 1].min()) / (p[:, 0].max() - p[:, 0].min()); RGB, EXT = terrain_rgb()
    W = 11.0; TOP, BOT = 1.6, 1.1; fig = plt.figure(figsize=(W, W * asp + TOP + BOT)); Hf = fig.get_figheight()
    ax = fig.add_axes([0.01, BOT / Hf, 0.98, W * asp / Hf], projection=PROJ); ax.set_extent(ext, crs=PC)
    draw_terrain(ax, RGB, EXT, regrid=2000); draw_base(ax)
    for o in offs:
        r = F[o]; im = imshow_proj(ax, r["f"], r["imext"], regrid_shape=2600, cmap=DIV, vmin=-2, vmax=2, interpolation="nearest", alpha=0.93, zorder=4)
        if r["XM"].any():
            imshow_proj(ax, np.where(r["XM"], 1.0, np.nan), r["imext"], regrid_shape=1600, cmap=ListedColormap(["#9a9a9a"]), vmin=0, vmax=1, alpha=0.55, interpolation="nearest", zorder=6)
            ax.contourf(r["CL"], r["CA"], r["XM"].astype(float), levels=[0.5, 1.5], colors="none", hatches=["////"], transform=PC, zorder=7)
    for o in offs:
        g = F[o]["geom"]; c = OFFICES[o]["color"]
        ax.add_geometries([g], PC, facecolor="none", edgecolor="white", linewidth=3.2, zorder=8)
        ax.add_geometries([g], PC, facecolor="none", edgecolor=c, linewidth=1.5, zorder=9)
        pt = cwa[o].representative_point()
        ax.text(pt.x, pt.y, f"{o}\n{F[o]['area']:+.2f} in", transform=PC, ha="center", va="center", fontsize=12, weight="bold", color="#111", zorder=12,
                path_effects=[pe.withStroke(linewidth=3.2, foreground="white")])
    xh = [Patch(fc="#9a9a9a", ec="#333", hatch="////", alpha=0.8, label="Excluded: radar beam-blockage artifact")] if any(F[o]["XM"].any() for o in offs) else []
    cartography(ax, ext, focus=None, shown=offs, extra=xh, reserve=[(cwa[o].representative_point().x, cwa[o].representative_point().y) for o in offs], states_over_data=True)
    vals = np.concatenate([F[o]["f"][np.isfinite(F[o]["f"])] for o in offs]); ext_ = {(0, 0): "neither", (1, 0): "min", (0, 1): "max", (1, 1): "both"}[(int(vals.min() < -2), int(vals.max() > 2))]
    cb = fig.colorbar(im, ax=ax, orientation="horizontal", fraction=0.04, pad=0.02, aspect=40, shrink=0.95, extend=ext_)
    cb.set_label("mean error, in per storm\nred = high bias (forecast too high);  blue = low bias (forecast too low)", fontsize=10.5); cb.ax.tick_params(labelsize=9.5)
    season = f", {int(season_id[:4])}–{str(int(season_id[:4])+1)[2:]} season" if season_id else ""
    d0 = min(F[o]["dates"][0] for o in offs); d1 = max(F[o]["dates"][1] for o in offs)
    fig.suptitle(f"New England NWS Offices{season}: Mean Snowfall Forecast Error", fontsize=19, weight="bold", y=1 - 0.12 / Hf, va="top")
    sub = "F = official NWS NDFD snowfall forecast issued ~12 h before onset;  O = NOHRSC snowfall analysis;  error = F − O.  Label = office area-wide mean error.\n" \
          "Storms: " + ",  ".join(f"{o} {F[o]['NS']}" for o in offs) + f"   ({d0:%b %Y} – {d1:%b %Y})"
    fig.text(0.5, 1 - 0.62 / Hf, sub, ha="center", va="top", fontsize=10.5, color="#444", linespacing=1.5)
    out = out_file or str(config.FIGS / (f"season_{season_id}/" if season_id else "") / "regional_mean_error_map.png")
    os.makedirs(os.path.dirname(out), exist_ok=True); fig.savefig(out, dpi=150, pad_inches=0.1); print("saved", out)

# ======================================================================================================================================
# RADAR BEAM-BLOCKAGE TEST (command: radar-test)
# ======================================================================================================================================
def radar_ray_test(off, lo0, la0, out_png=None):
    """Radar beam-blockage check.  Flank-relative azimuthal anomaly of ln(O/F) around a radar site.  A narrow azimuth band whose cells are all
    low (or high) relative to neighbouring azimuths at the same range, for >10 consecutive range gates, is a radar-geometry artifact in the
    NOHRSC analysis.  Use the printed azimuth profile to choose an entry for RADAR_EXCLUSIONS."""
    D = pickle.load(open(config.STORMS / f"ndfd_{off}.pkl", "rb")); R = D["R"]; KEYS = sorted(R); CL, CA, CM = D["CLON"], D["CLAT"], D["CMASK"]
    F=np.stack([R[k]["fc"] for k in KEYS]); T=np.stack([R[k]["nsa"] for k in KEYS]); E=F-T
    ok=np.isfinite(E)&CM[None]
    with np.errstate(all="ignore"):
        mF=np.nanmean(np.where(ok,F,np.nan),0); mO=np.nanmean(np.where(ok,T,np.nan),0); mE=mF-mO
    LON,LAT=np.meshgrid(CL,CA)
    def anomaly(lo0,la0):
        dy=(LAT-la0)*111.0; dx=(LON-lo0)*111.0*np.cos(np.radians(la0)); rng=np.hypot(dx,dy); az=(np.degrees(np.arctan2(dx,dy))+360)%360
        V=CM&np.isfinite(mE)&(mF>0.3); lr=np.where(V,np.log((mO+0.3)/(mF+0.3)),np.nan)
        ii=np.argwhere(V); r_=rng[V]; a_=az[V]; l_=lr[V]; an=np.full(len(ii),np.nan)
        for n,(r,a) in enumerate(zip(r_,a_)):
            if r<15 or r>130: continue
            da=np.abs((a_-a+180)%360-180); arc=da*np.pi/180*r
            m=(np.abs(r_-r)<=3.0)&(arc>=6)&(arc<=24)
            if m.sum()>=6: an[n]=l_[n]-np.median(l_[m])
        A=np.full(CM.shape,np.nan); A[tuple(ii.T)]=an; return A,rng,az
    A,rng,az=anomaly(lo0,la0)
    # ray persistence: for each 3-deg azimuth bin, fraction of valid cells with anomaly < -0.12 (high bias dip) and > +0.12, ranges 20-120 km
    for sgn,name in ((-1,"dip (O lower than flanks)"),(1,"bump (O higher than flanks)")):
        out=[]
        for a0 in range(0,360,3):
            m=np.isfinite(A)&(az>=a0)&(az<a0+3)&(rng>20)&(rng<120)
            if m.sum()>=8: out.append((float(np.mean(sgn*A[m]>0.12)),a0,int(m.sum()),float(rng[m].max())))
        out.sort(reverse=True); print(name,"worst azimuth bins (frac cells beyond 0.12, az, n, max range km):",[(round(f,2),a,n,round(r)) for f,a,n,r in out[:6]])
    fig,ax=plt.subplots(1,2,figsize=(18,7)); ex=[CL[0]-.02,CL[-1]+.02,CA[-1]-.02,CA[0]+.02]; asp=1/np.cos(np.radians(la0))
    ax[0].imshow(np.where(CM,mE,np.nan),extent=ex,cmap="RdBu_r",vmin=-2,vmax=2,aspect=asp); ax[1].imshow(A,extent=ex,cmap="PuOr",vmin=-0.5,vmax=0.5,aspect=asp)
    for a_ in ax:
        a_.plot(lo0,la0,"k^"); 
        for azd in range(0,360,30):
            a_.plot([lo0,lo0+120*np.sin(np.radians(azd))/(111*np.cos(np.radians(la0)))],[la0,la0+120*np.cos(np.radians(azd))/111],color="k",lw=0.3,alpha=0.4); 
        a_.set_xlim(ex[0],ex[1]); a_.set_ylim(ex[2],ex[3])
    ax[0].set_title("mean error F−O (red = high bias)"); ax[1].set_title("azimuthal anomaly ln(O/F) vs flanks (orange = O lower)")
    fig.savefig(out_png or str(config.FIGS/"radar_ray_test.png"),dpi=65,bbox_inches="tight")
    # 1-degree azimuth profile around the worst dip, to pick the exclusion sector edges
    worst=None
    for a0 in range(0,360):
        m=np.isfinite(A)&(az>=a0)&(az<a0+1)&(rng>20)&(rng<110)
        if m.sum()>=4:
            v=float(np.median(A[m]))
            if worst is None or v<worst[0]: worst=(v,a0)
    print("most negative 1-deg azimuth bin (median anomaly, azimuth):",worst)
    if worst:
        row=[]
        for a0 in range(worst[1]-12,worst[1]+13):
            a=a0%360; m=np.isfinite(A)&(az>=a)&(az<a+1)&(rng>20)&(rng<110)
            row.append((a,round(float(np.median(A[m])),2) if m.sum()>=4 else None))
        print("profile around it (azimuth, median anomaly):",row)


# ======================================================================================================================================
# COMMAND LINE / ORCHESTRATION
# ======================================================================================================================================
ALL_STAGES = ["terrain", "cube", "storms", "figures"]
FIG_STAGES = {  # code: (output file suffix, description)
    "e": ("_0_mean_error_map.png",           "large Mean Error map"),
    "a": ("_1_bias_MAE_maps.png",            "Mean Error, MAE, Percent Bias, Normalized MAE"),
    "b": ("_2_hit_rate_maps.png",            "Hit Rate, Forecast Tendency, Reliability Category"),
    "c": ("_3_error_charts.png",             "hit/high/low rates by amount and terrain; F vs O"),
    "d": ("_4_bias_terrain_coast.png",       "Percent Bias vs elevation and distance to ocean"),
    "g": ("_5_largest_storms_F_O_error.png", "four largest storms: F, O, error"),
    "m": ("_6_mean_F_O_error_maps.png",      "mean forecast, mean observed, mean error"),
}

def preflight(stages):
    need = {"cube": ["tifffile"], "storms": ["pygrib"]}; missing = []
    for s in stages:
        for mod in need.get(s, []):
            try: __import__(mod)
            except Exception as e: missing.append(f"{mod} ({e.__class__.__name__}) - needed by the '{s}' stage")
    if missing: sys.exit("Missing Python packages: " + "; ".join(missing) + "\nInstall with: pip install tifffile pygrib   (pygrib needs the ECCODES library)")
    free = shutil.disk_usage(config.WORK if config.WORK.is_dir() else ".").free / 1e9
    if free < 12: print(f"WARNING: only {free:.0f} GB free; the full record needs about 8 GB of working space.", flush=True)

def default_jobs():
    """parallel processes: CPU count, limited by available RAM (a figure process peaks at ~2.5 GB)"""
    mem = 8.0
    try: mem = int([l for l in open("/proc/meminfo") if l.startswith("MemAvailable")][0].split()[1]) / 1e6
    except Exception: pass
    return max(1, min(os.cpu_count() or 1, int(mem // 3), 4))

def banner(s): print(f"\n=== {s}  [{dt.datetime.now():%H:%M:%S}]", flush=True)

def run(cmd, env=None, log=None):
    t0 = time.time()
    with open(log, "a") if log else open(os.devnull, "w") as lf:
        r = subprocess.run(cmd, env=env, stdout=lf, stderr=subprocess.STDOUT)
    return r.returncode, time.time() - t0

def _worker(argv):
    """internal entry points: the orchestrator runs one figure process per office (memory isolation, shared caches) and one for the regional map"""
    kind = argv[0]; ap = argparse.ArgumentParser(prog=f"{SELF.name} {kind}")
    if kind == "_figures":
        ap.add_argument("office"); ap.add_argument("stages"); ap.add_argument("--outdir", required=True); ap.add_argument("--season", default="")
        a = ap.parse_args(argv[1:]); sys.exit(0 if draw_office_figures(a.office, a.stages.split(","), a.outdir, a.season) else 3)
    if kind == "_regional":
        ap.add_argument("--offices", nargs="+"); ap.add_argument("--season", default=None); ap.add_argument("--out", default=None)
        a = ap.parse_args(argv[1:]); regional_map(a.offices, a.season, a.out); sys.exit(0)
    if kind == "radar-test":
        ap.add_argument("office"); ap.add_argument("radar_lon", type=float); ap.add_argument("radar_lat", type=float); ap.add_argument("out", nargs="?")
        a = ap.parse_args(argv[1:]); radar_ray_test(a.office.upper(), a.radar_lon, a.radar_lat, a.out); sys.exit(0)

def main():
    if len(sys.argv) > 1 and sys.argv[1] in ("_figures", "_regional", "radar-test"): return _worker(sys.argv[1:])
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--workdir", default="work", help="directory for all downloads and outputs (default ./work)")
    ap.add_argument("--offices", nargs="+", default=None, help="subset of the configured offices (default all)")
    ap.add_argument("--stages", default=",".join(ALL_STAGES), help="comma list from: " + ",".join(ALL_STAGES))
    ap.add_argument("--figures", default=None, help="figure types to draw; default all (e,a,m,b,g,c,d), or e,a,m with --season")
    ap.add_argument("--season", default=None, help="verify one winter only, e.g. 2025-26 (storms starting Jul 2025 - Jun 2026); figures go to figures/season_2025-26/")
    ap.add_argument("--start", default=None, help=f"first date, YYYY-MM-DD (default {DEFAULT_START}; the NDFD archive begins June 2018)")
    ap.add_argument("--end", default=None, help=f"last date, YYYY-MM-DD (default {DEFAULT_END})")
    ap.add_argument("--workers", type=int, default=6, help="parallel download threads (the per-host connection/rate caps in NET_LIMITS still apply)")
    ap.add_argument("--jobs", type=int, default=None, help="parallel processes for decoding and figure drawing (default: CPUs, limited by RAM)")
    ap.add_argument("--nohrsc-rate", type=float, default=None, help="max requests/second to nohrsc.noaa.gov (default 8)")
    ap.add_argument("--ncei-rate", type=float, default=None, help="max requests/second to ncei.noaa.gov (default 3)")
    ap.add_argument("--cwa-file", default=None, help="GeoJSON with the CWA polygons (property 'cwa'); default: downloaded once from the NWS reference map service")
    ap.add_argument("--force", action="store_true", help="redo stages even if their outputs exist")
    a = ap.parse_args()

    os.environ["NWSNOW_WORKDIR"] = str(pathlib.Path(a.workdir).resolve())
    if a.start: os.environ["NWSNOW_START"] = a.start
    if a.end: os.environ["NWSNOW_END"] = a.end
    if a.nohrsc_rate: os.environ["NWSNOW_NOHRSC_RATE"] = str(a.nohrsc_rate)
    if a.ncei_rate: os.environ["NWSNOW_NCEI_RATE"] = str(a.ncei_rate)
    if a.cwa_file: os.environ["NWSNOW_CWA_FILE"] = str(pathlib.Path(a.cwa_file).resolve())
    _init_config(); config.WORK.mkdir(parents=True, exist_ok=True)
    offices = [o.upper() for o in (a.offices or config.ORDER)]
    bad = [o for o in offices if o not in OFFICES]
    if bad: sys.exit(f"Unknown office(s): {bad}.  Choose from {list(OFFICES)}")
    stages = [s.strip() for s in a.stages.split(",") if s.strip()]
    for s in stages:
        if s not in ALL_STAGES: sys.exit(f"Unknown stage '{s}'. Choose from {ALL_STAGES}")
    figs = [f.strip() for f in (a.figures or ("e,a,m" if a.season else "e,a,m,b,g,c,d")).split(",") if f.strip()]
    for f in figs:
        if f not in FIG_STAGES: sys.exit(f"Unknown figure type '{f}'. Choose from {list(FIG_STAGES)}")
    figdir = config.FIGS / f"season_{a.season}" if a.season else config.FIGS
    preflight(stages)
    for d in (config.LOGS, figdir): d.mkdir(parents=True, exist_ok=True)
    if "storms" in stages or "figures" in stages: ensure_cwa_file()
    print(f"work directory : {config.WORK}\nperiod         : {config.START:%Y-%m-%d} .. {config.END:%Y-%m-%d}\noffices        : {offices}\nstages         : {stages}", flush=True)
    jobs = a.jobs or default_jobs(); print(f"parallel jobs  : {jobs}", flush=True)
    T0 = time.time()

    if "terrain" in stages:
        banner("terrain: AWS terrain tiles -> dem.npz")
        if config.DEM.exists() and not a.force: print("  have dem.npz, skipping", flush=True)
        else: stage_terrain()
    if "cube" in stages:
        banner("cube: NOHRSC 6-h snowfall (New England window)")
        stage_cube(workers=a.workers)
    if "storms" in stages:
        banner("storms: detect storms, fetch NDFD forecasts")
        todo = [o for o in offices if a.force or not (config.STORMS / f"ndfd_{o}.pkl").exists()]
        if not todo: print("  all storm files exist, skipping (use --force to rebuild)", flush=True)
        else: stage_storms(todo, workers=a.workers, jobs=jobs)
    if "figures" in stages:
        banner("figures")
        def draw(o):
            todo = [f for f in figs if a.force or not (figdir / f"{o}{FIG_STAGES[f][0]}").exists()]
            if not todo: return f"  {o}: all {len(figs)} figures exist"
            cmd = [sys.executable, str(SELF), "_figures", o, ",".join(todo), "--outdir", str(figdir)] + (["--season", a.season] if a.season else [])
            rc, secs = run(cmd, env=dict(os.environ), log=config.LOGS / f"figures_{o}.log")      # one process per office draws all its figures
            what = "ok" if rc == 0 else ("no storms in this season" if rc == 3 else f"FAILED (see logs/figures_{o}.log)")
            return f"  {o}: {len(todo)} figures {what}  ({secs:.0f} s)"
        ready = [o for o in offices if (config.STORMS / f"ndfd_{o}.pkl").exists()]
        for o in offices:
            if o not in ready: print(f"  {o}: no storm file - run the 'storms' stage first", flush=True)
        if ready:
            config.TERRAIN.mkdir(parents=True, exist_ok=True)
            first = draw(ready[0]) if not ((config.TERRAIN / "basemap_lcc_v2.pkl").exists() and (config.TERRAIN / "reference_layers.pkl").exists()) else None   # builds the cached base/reference layers once, before going parallel
            if first: print(first, flush=True)
            with ThreadPoolExecutor(jobs) as ex:
                for msg in ex.map(draw, ready[1 if first else 0:]): print(msg, flush=True)
            reg = [o for o in ready if o not in config.REGIONAL_EXCLUDE]       # one map with all offices on a common color scale
            regout = figdir / "regional_mean_error_map.png"
            if len(reg) >= 2 and (a.force or not regout.exists()):
                cmd = [sys.executable, str(SELF), "_regional", "--offices", *reg, "--out", str(regout)] + (["--season", a.season] if a.season else [])
                rc, secs = run(cmd, env=dict(os.environ), log=config.LOGS / "regional.log")
                print(f"  regional map ({' '.join(reg)}): {'ok' if rc == 0 else 'FAILED (see logs/regional.log)'}  ({secs:.0f} s)", flush=True)
    print(f"\nfinished in {(time.time() - T0) / 60:.1f} min.  Figures: {figdir}", flush=True)

if __name__ == "__main__":
    main()
