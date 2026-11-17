# Snowfall forecast error viewer

A small web page that shows the **mean snowfall forecast error (forecast − observed) for each New England NWS office** on an interactive map. The error grids are plain **`.bin` files** read by the browser; PHP serves them, JavaScript paints them.

## Run it

```bash
# 1. write the grids for the most recent season (needs the storm files from the pipeline: work/storms/ndfd_<OFFICE>.pkl)
python scripts/export_viewer_data.py --workdir ./work                 # --season 2025-26 to pick one explicitly

# 2. serve the viewer (PHP 8.1 or newer)
php -S localhost:8000 -t viewer                                       # then open http://localhost:8000
```

The repository already contains the exported 2025–26 season (`viewer/data/`, about 160 KB), so step 2 alone works after cloning. Running step 1 again for another season adds it to the season menu.

### GitHub Pages (no PHP)

GitHub Pages cannot run PHP, so the repository also contains a static twin of the page, `viewer/index.html`, which reads `data/manifest.json`, `data/cwa.geojson` and the `.bin` files directly. The workflow `.github/workflows/pages.yml` publishes it on every push to `main` that touches `viewer/`.

One-time setup: **Settings → Pages → Build and deployment → Source: GitHub Actions**. The site then appears at `https://<user>.github.io/<repo>/` (the first run can also be started from the Actions tab). To try the static version locally: `python -m http.server -d viewer` and open `/index.html`.

`index.php` and `index.html` are the same page; if you edit one, edit the other. `app.js` is shared and switches mode on `window.NWS.dataDir`.

### Any PHP host

On a normal web host, upload the `viewer/` folder to any PHP 8.1+ site (Apache or nginx + PHP-FPM). `.htaccess` turns off directory listing. To keep the data outside the web root, set the environment variable `NWS_VIEWER_DATA=/path/to/data`.

## What it shows

* **Overlay:** the mean error of each ~4 km (0.04°) cell over the season's storms, drawn **unsmoothed at the native grid** on the same fixed diverging scale as the printed maps (±2 in; red = forecast too high, blue = too low). Cells are painted by sampling the grid per screen pixel, so they stay square and correctly placed at every zoom level.
* **CWA outlines** (colored per office) with the office's area-wide mean error as a label. Click an outline, or a row in the table, for details.
* **Details:** storms, mean error with its 90 % storm-bootstrap interval, percent bias, MAE, RMSE, hit rate, share of the CWA above +0.5 in / below −0.5 in, and the full-record mean error for comparison. Offices with fewer than 10 storms are flagged.
* **Cursor tooltip:** a small window follows the pointer over the map and reads the office and mean error of the cell under it (tap on touch screens).
* **Hatched cells** are radar beam-blockage artifacts in the observation analysis, excluded from every statistic.
* **Own basemap** (see below): shaded-relief terrain, ocean and lakes, coastline, the U.S.–Canada border, state and province lines, and place names. Nothing is loaded from third-party tile servers.
* Controls: season, basemap (terrain or plain), overlay opacity, show/hide per office. The state is kept in the URL hash, so a view can be shared: `#season=2025-26&office=CAR&base=terrain&op=0.7`. `Esc` clears the selection.

## Basemap (`viewer/basemap/`)

Built once by `python scripts/build_basemap.py` (about 40 s; needs internet only for that step) and committed, so the viewer itself is self-contained.

| File | Content |
|---|---|
| `tiles/{z}/{x}/{y}.webp` | shaded-relief terrain, Web Mercator XYZ tiles, zoom 5–9 (the map allows zoom 10 by stretching zoom-9 tiles) |
| `water.geojson` | ocean and lake polygons |
| `lines.geojson` | coastline and lake shores, the international border, state/province lines |
| `labels.json` | state, lake and ocean names with the zoom they appear at |
| `meta.json` | extent (39.3–48.7°N, 77.5–65.0°W), zoom range, sources |

* **Terrain:** elevation from *Terrain Tiles on AWS* (Terrarium PNGs on S3, `elevation-tiles-prod`), hillshaded by the script (light from the northwest, with the ground size of a pixel corrected for latitude) and tinted by elevation in muted colours. Data sources of those tiles (USGS NED, NRCan CDEM, SRTM and others) are listed at <https://github.com/tilezen/joerd/blob/master/docs/attribution.md>.
* **Boundaries, coast, lakes:** Natural Earth 1:50m (medium resolution; public domain), clipped to the extent and rounded to 0.0001°. They are vectors, so lines stay sharp at every zoom. Draw order: terrain, water, snowfall-error overlay, coast / state / border lines, CWA outlines, place names, office labels.
* Options: `--zmax 10` adds a zoom level (about four times the tiles), `--west/--east/--south/--north` change the extent, `--quality` the WebP quality.

## Data files (`viewer/data/`)

| File | Content |
|---|---|
| `manifest.json` | what exists: layer definition (color range), offices and colors, seasons, and for each office the file name, grid geometry and statistics |
| `<season>/<OFFICE>_mean_error.bin` | the grid |
| `cwa.geojson` | simplified CWA outlines (EPSG:4326) |

**`.bin` format.** Raw **float32, little-endian, row-major, no header**. Row 0 is the **north** edge. Geometry comes from the manifest (`nx`, `ny`, `west`, `north`, `dx`, `dy`, in degrees; `west` and `north` are the outer edges of the first column and row). Special values: **`NaN`** = outside the CWA or fewer than 90 % of storms have data; **`-999`** = excluded radar sector.

```python
import json, numpy as np
m = json.load(open("viewer/data/manifest.json")); e = m["seasons"]["2025-26"]["offices"]["CAR"]
grid = np.fromfile("viewer/data/" + e["file"], dtype="<f4").reshape(e["ny"], e["nx"])      # inches per storm
lon = e["west"]  + (np.arange(e["nx"]) + 0.5) * e["dx"]                                        # cell-center longitudes
lat = e["north"] - (np.arange(e["ny"]) + 0.5) * e["dy"]                                        # cell-center latitudes (north to south)
```

The values are exactly those of the printed mean-error maps (pipeline figure stage `e`) and the statistics match `results/office_summary.csv`.

## Files

```
viewer/
  index.php          page shell for PHP hosting; inlines the manifest
  index.html         the same page for static hosting (GitHub Pages)
  api.php            read-only API: ?action=manifest | cwa | bin&season=2025-26&office=CAR
  config.php         title and data directory
  assets/app.js      map, basemap layers, .bin decoding, tile painting, table, legend, cursor tooltip, details
  assets/style.css
  vendor/leaflet/    Leaflet 1.9.4 (BSD-2-Clause), so the viewer needs no CDN for its code
  data/              exported grids, manifest, outlines
  basemap/           own basemap: terrain tiles, water, lines, labels (scripts/build_basemap.py)
```

`api.php` looks the season and office up as **keys of the manifest** and never builds a file path from request input, so only the listed files can be served. Responses carry an `ETag` (304 when unchanged), are gzip-compressed when the client accepts it, and are cacheable for an hour. Only `GET` and `HEAD` are allowed.

## Notes

* The viewer makes **no requests to other servers**: the basemap, the overlay and the JavaScript library are all served from this folder (an earlier version used CARTO tiles, which now require an API key).
* A single season rests on few storms (5–24 per office), so differences between offices in one season are often within the intervals shown. See the main README for the multi-season study.
* Data: NWS NDFD forecasts (NCEI), NOHRSC National Snowfall Analysis, NWS county warning areas. Leaflet © Vladimir Agafonkin and contributors.
