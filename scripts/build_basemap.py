#!/usr/bin/env python3
"""Build the web viewer's own basemap, so the viewer needs no third-party tile server.

    python scripts/build_basemap.py                        # writes viewer/basemap/ (about 20-30 s after the first download)

What it makes (viewer/basemap/):
    tiles/{z}/{x}/{y}.webp   shaded-relief terrain tiles (standard XYZ / Web Mercator, 256 px), zoom 5-9
    water.geojson            ocean and lake polygons                 (Natural Earth 50m)
    lines.geojson            coastline + lake shores, the U.S.-Canada border, state/province lines  (Natural Earth 50m)
    labels.json              state/province, lake and ocean name labels with a minimum zoom
    meta.json                extent, zoom range, tile size, sources

Terrain: "Terrain Tiles on AWS" (Terrarium-encoded elevation PNGs, s3://elevation-tiles-prod), hillshaded here from the elevation itself
(light from the northwest, cell size corrected for the Web Mercator scale at each latitude) and tinted by elevation in muted colours so the
snowfall-error overlay stays the focus.  Boundaries: Natural Earth 1:50m (medium resolution, public domain), clipped to the extent and
rounded to 0.0001 degree.  Downloads are cached in --cache, so re-running is quick.
"""
import argparse, io, json, math, sys, time, urllib.request, warnings
from concurrent.futures import ThreadPoolExecutor
from datetime import date
from pathlib import Path

import numpy as np
from PIL import Image

TERRARIUM = "https://s3.amazonaws.com/elevation-tiles-prod/terrarium/{z}/{x}/{y}.png"
UA = "ne-wfo-snowfall-viewer-basemap/1.0 (+https://github.com/aaronperry-dev/NWS-WFO-Snowfall-Error-Viewer)"
R_EARTH = 6378137.0

# elevation (m) -> land colour; deliberately pale and low-chroma
TINT = [(0, (205, 220, 189)), (150, (219, 228, 194)), (400, (232, 226, 193)), (800, (222, 208, 178)), (1300, (210, 198, 180)), (2000, (232, 228, 222))]
# state-name positions moved off the default label point where it would sit under the viewer's office labels (BOX, OKX)
LABEL_AT = {"Massachusetts": (-72.55, 42.38), "Connecticut": (-72.75, 41.62), "New York": (-75.45, 42.75)}


def tx(lon, z): return (lon + 180.0) / 360.0 * 2 ** z
def ty(lat, z): return (1.0 - math.asinh(math.tan(math.radians(lat))) / math.pi) / 2.0 * 2 ** z
def lat_of_ty(t, z): return math.degrees(math.atan(math.sinh(math.pi * (1.0 - 2.0 * t / 2 ** z))))


# ---------------------------------------------------------------- terrain
def fetch_tile(cache, z, x, y):
    n = 2 ** z; x %= n; y = min(max(y, 0), n - 1)
    f = cache / "terrarium" / str(z) / str(x) / f"{y}.png"
    if f.exists() and f.stat().st_size > 0: return f
    f.parent.mkdir(parents=True, exist_ok=True)
    for attempt in range(4):
        try:
            req = urllib.request.Request(TERRARIUM.format(z=z, x=x, y=y), headers={"User-Agent": UA})
            with urllib.request.urlopen(req, timeout=60) as r: data = r.read()
            f.write_bytes(data); return f
        except Exception as e:
            if attempt == 3: raise RuntimeError(f"tile {z}/{x}/{y}: {e}")
            time.sleep(1.5 * (attempt + 1))


def elevation_mosaic(cache, z, x0, x1, y0, y1, pool):
    """elevation (m) of tiles x0-1..x1+1, y0-1..y1+1 (one tile of margin so edges shade correctly)"""
    xs, ys = range(x0 - 1, x1 + 2), range(y0 - 1, y1 + 2)
    files = list(pool.map(lambda xy: fetch_tile(cache, z, *xy), [(x, y) for y in ys for x in xs]))
    m = np.zeros((len(ys) * 256, len(xs) * 256), np.float32); k = 0
    for j in range(len(ys)):
        for i in range(len(xs)):
            a = np.asarray(Image.open(files[k]).convert("RGB"), np.float32); k += 1
            m[j * 256:(j + 1) * 256, i * 256:(i + 1) * 256] = a[..., 0] * 256 + a[..., 1] + a[..., 2] / 256 - 32768
    return m


def shade_mosaic(elev, z, y_first, exag=2.0, azimuth=315.0, altitude=42.0, strength=1.25):
    """RGB uint8 image: elevation tint x hillshade"""
    h, w = elev.shape; zz = np.maximum(elev, 0.0)                         # no bathymetry (the viewer paints the water on top)
    rows = (np.arange(h) + 0.5) / 256.0 + y_first                          # fractional tile-y of each pixel row
    lat = np.array([lat_of_ty(t, z) for t in rows])
    cell = (2 * math.pi * R_EARTH / (256.0 * 2 ** z) * np.cos(np.radians(lat)))[:, None]   # metres per pixel on the ground
    p = np.gradient(zz, axis=1) / cell * exag                               # dz/dx (east)
    q = -np.gradient(zz, axis=0) / cell * exag                              # dz/dy (north; image rows run south)
    az, alt = math.radians(azimuth), math.radians(altitude); zen = math.pi / 2 - alt
    lx, ly, lz = math.sin(zen) * math.sin(az), math.sin(zen) * math.cos(az), math.cos(zen)
    illum = (-p * lx - q * ly + lz) / np.sqrt(1 + p * p + q * q)
    factor = np.clip(1 + strength * (illum / lz - 1), 0.42, 1.3)           # 1.0 on flat ground
    el = np.array([t[0] for t in TINT], float)
    base = np.stack([np.interp(zz, el, [t[1][c] for t in TINT]) for c in range(3)], -1)
    return np.clip(base * factor[..., None], 0, 255).astype(np.uint8)


def build_tiles(args, out, cache, pool):
    W, E, S, N = args.west, args.east, args.south, args.north
    total = 0; size = 0
    for z in range(args.zmin, args.zmax + 1):
        x0, x1, y0, y1 = int(tx(W, z)), int(tx(E, z)), int(ty(N, z)), int(ty(S, z))
        t0 = time.time(); elev = elevation_mosaic(cache, z, x0, x1, y0, y1, pool)
        rgb = shade_mosaic(elev, z, y0 - 1)
        n = 0
        for y in range(y0, y1 + 1):
            for x in range(x0, x1 + 1):
                j, i = (y - y0 + 1) * 256, (x - x0 + 1) * 256
                f = out / "tiles" / str(z) / str(x) / f"{y}.webp"; f.parent.mkdir(parents=True, exist_ok=True)
                Image.fromarray(rgb[j:j + 256, i:i + 256]).save(f, "WEBP", quality=args.quality, method=6); size += f.stat().st_size; n += 1
        total += n; print(f"  z{z}: {n:4d} tiles  ({time.time() - t0:.0f} s)")
    return total, size


# ---------------------------------------------------------------- vectors
def build_vectors(args, out, cache):
    import cartopy, shapefile
    import cartopy.io.shapereader as shp
    from shapely import set_precision
    from shapely.geometry import box, mapping, shape
    cartopy.config["data_dir"] = str(cache / "cartopy")
    warnings.filterwarnings("ignore")
    W, E, S, N = args.west, args.east, args.south, args.north
    ext, clip = box(W, S, E, N), box(W - 0.6, S - 0.6, E + 0.6, N + 0.6)

    def records(cat, name):
        r = shapefile.Reader(shp.natural_earth("50m", cat, name))
        keys = [f[0].lower() for f in r.fields[1:]]
        for sr in r.shapeRecords():
            if sr.shape.shapeType == 0: continue
            yield shape(sr.shape.__geo_interface__), dict(zip(keys, sr.record))

    def feat(g, **props):
        g = set_precision(g, 1e-4)
        return {"type": "Feature", "properties": props, "geometry": mapping(g)} if not g.is_empty else None

    def lines_in(g):
        from shapely.geometry import LineString, MultiLineString
        g = g.intersection(clip)
        parts = [g] if g.geom_type == "LineString" else [p for p in getattr(g, "geoms", []) if p.geom_type == "LineString"]
        return [p for p in parts if p.length > 0.01]

    # water polygons
    water = []
    for g, r in records("physical", "ocean"):
        if g.intersects(clip): water.append(feat(g.intersection(clip), kind="ocean"))
    lakes = [(g, r) for g, r in records("physical", "lakes") if g.intersects(clip)]
    for g, r in lakes: water.append(feat(g.intersection(clip), kind="lake", name=r.get("name") or ""))
    water = [f for f in water if f]

    # lines: coast (ocean + lake shores), international border, state/province lines
    lines = []
    for g, r in records("physical", "coastline"):
        if g.intersects(clip): lines += [feat(p, kind="coast") for p in lines_in(g)]
    for g, r in lakes: lines += [feat(p, kind="coast") for p in lines_in(g.boundary)]
    border = None
    for g, r in records("cultural", "admin_0_boundary_lines_land"):
        if g.intersects(clip):
            lines += [feat(p, kind="border") for p in lines_in(g)]
            border = g if border is None else border.union(g)
    mask = border.buffer(0.02) if border is not None else None
    for g, r in records("cultural", "admin_1_states_provinces_lines"):
        if not g.intersects(clip): continue
        if mask is not None: g = g.difference(mask)                          # the country border is drawn once, as the border
        if not g.is_empty: lines += [feat(p, kind="state") for p in lines_in(g)]
    lines = [f for f in lines if f]

    # labels
    labels = []
    for g, r in records("cultural", "admin_1_states_provinces"):
        if not g.intersects(ext): continue
        part = g.intersection(ext)
        if part.area < 0.25: continue
        pt = max(part.geoms, key=lambda q: q.area) if part.geom_type == "MultiPolygon" else part
        c = pt.representative_point(); mz = 5 if part.area > 8 else 6 if part.area > 2.5 else 7 if part.area > 0.9 else 8
        lon, lat = LABEL_AT.get(r.get("name"), (c.x, c.y))
        labels.append({"name": r.get("name"), "lon": round(lon, 3), "lat": round(lat, 3), "kind": "state", "minzoom": mz})
    for g, r in lakes:
        part = g.intersection(ext); nm = r.get("name")
        if not nm or part.is_empty or part.area < 0.015 or "servoir" in nm or "eservoir" in nm: continue      # skip remote reservoirs
        pt = max(part.geoms, key=lambda q: q.area) if part.geom_type == "MultiPolygon" else part
        c = pt.representative_point(); mz = 6 if part.area > 0.4 else 8 if part.area > 0.08 else 9
        labels.append({"name": nm, "lon": round(c.x, 3), "lat": round(c.y, 3), "kind": "water", "minzoom": mz})
    labels += [{"name": "Atlantic Ocean", "lon": -66.6, "lat": 41.4, "kind": "ocean", "minzoom": 5},     # fixed positions: the ocean polygon has no useful label point
               {"name": "Gulf of Maine", "lon": -68.3, "lat": 43.1, "kind": "ocean", "minzoom": 7}]

    dump = lambda name, obj: (out / name).write_text(json.dumps(obj, separators=(",", ":"), ensure_ascii=False))
    dump("water.geojson", {"type": "FeatureCollection", "features": water})
    dump("lines.geojson", {"type": "FeatureCollection", "features": lines})
    dump("labels.json", labels)
    kinds = {k: sum(1 for f in lines if f["properties"]["kind"] == k) for k in ("coast", "border", "state")}
    print(f"  water: {len(water)} polygons ({len(lakes)} lakes) | lines: {kinds} | labels: {len(labels)}")


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--out", default="viewer/basemap"); ap.add_argument("--cache", default="work/basemap", help="download cache")
    ap.add_argument("--west", type=float, default=-77.5); ap.add_argument("--east", type=float, default=-65.0)
    ap.add_argument("--south", type=float, default=39.3); ap.add_argument("--north", type=float, default=48.7)
    ap.add_argument("--zmin", type=int, default=5); ap.add_argument("--zmax", type=int, default=9)
    ap.add_argument("--quality", type=int, default=80, help="WebP quality of the terrain tiles")
    ap.add_argument("--threads", type=int, default=6, help="parallel downloads (kept modest)")
    args = ap.parse_args()
    out, cache = Path(args.out), Path(args.cache); out.mkdir(parents=True, exist_ok=True); cache.mkdir(parents=True, exist_ok=True)
    print("terrain tiles ...")
    with ThreadPoolExecutor(args.threads) as pool: n, size = build_tiles(args, out, cache, pool)
    print("boundaries ...")
    build_vectors(args, out, cache)
    meta = {"bounds": [[args.south, args.west], [args.north, args.east]], "minzoom": args.zmin, "maxzoom": args.zmax, "tile_size": 256,
            "tiles": "tiles/{z}/{x}/{y}.webp", "built": date.today().isoformat(),
            "sources": {"terrain": "Terrain Tiles on AWS (Terrarium), hillshaded locally", "boundaries": "Natural Earth 1:50m"}}
    (out / "meta.json").write_text(json.dumps(meta, indent=1))
    print(f"done: {n} tiles, {size / 1e6:.1f} MB  ->  {out}")


if __name__ == "__main__":
    sys.exit(main())
