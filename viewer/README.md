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
* **Hover** (or tap on touch screens) reads the value of the cell under the pointer.
* **Hatched cells** are radar beam-blockage artifacts in the observation analysis, excluded from every statistic.
* Controls: season, basemap (CARTO light, OpenStreetMap, Esri topographic, Esri satellite), overlay opacity, show/hide per office. The state is kept in the URL hash, so a view can be shared: `#season=2025-26&office=CAR&base=topo&op=0.7`. `Esc` clears the selection.

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
  assets/app.js      map, .bin decoding, tile painting, table, legend, hover, details
  assets/style.css
  vendor/leaflet/    Leaflet 1.9.4 (BSD-2-Clause), so the viewer needs no CDN for its code
  data/              exported grids, manifest, outlines
```

`api.php` looks the season and office up as **keys of the manifest** and never builds a file path from request input, so only the listed files can be served. Responses carry an `ETag` (304 when unchanged), are gzip-compressed when the client accepts it, and are cacheable for an hour. Only `GET` and `HEAD` are allowed.

## Notes

* **Basemaps** are loaded by the visitor's browser from CARTO, OpenStreetMap and Esri, subject to their terms and usage policies; the OpenStreetMap standard tile server is meant for light use. For a busy public site use a provider with an account and key. The overlay itself needs no external service.
* A single season rests on few storms (5–24 per office), so differences between offices in one season are often within the intervals shown. See the main README for the multi-season study.
* Data: NWS NDFD forecasts (NCEI), NOHRSC National Snowfall Analysis, NWS county warning areas. Leaflet © Vladimir Agafonkin and contributors.
