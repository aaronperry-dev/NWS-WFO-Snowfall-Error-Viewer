/* NWS WFO snowfall forecast error viewer.
 * Loads one raw float32 .bin grid per office (via api.php, or straight from data/ on static hosting), and paints them on a Leaflet map with a custom tile layer that samples the
 * grids per screen pixel (nearest cell), so the cells stay crisp, unsmoothed and correctly placed at every zoom.
 * Grid file format (see manifest.json): little-endian float32, row-major, row 0 = north edge, NaN = no data, -999 = excluded radar sector. */
(() => {
'use strict';
/* Two ways to get the data: through api.php (PHP hosting), or, when the page sets window.NWS.dataDir, straight from static files (GitHub Pages). */
const API = (window.NWS && window.NWS.api) || 'api.php', DATA_DIR = window.NWS && window.NWS.dataDir;
const url = {
  manifest: () => DATA_DIR ? `${DATA_DIR}manifest.json` : `${API}?action=manifest`,
  cwa: () => DATA_DIR ? `${DATA_DIR}cwa.geojson` : `${API}?action=cwa`,
  bin: (season, off, e) => DATA_DIR ? DATA_DIR + e.file.split('/').map(encodeURIComponent).join('/') : `${API}?action=bin&season=${encodeURIComponent(season)}&office=${encodeURIComponent(off)}`,
};
const MINUS = '−', EXCLUDED = -999;
const el = (id) => document.getElementById(id);
const esc = (s) => String(s).replace(/[&<>"']/g, (c) => ({ '&': '&amp;', '<': '&lt;', '>': '&gt;', '"': '&quot;', "'": '&#39;' }[c]));
const fmt = (v, d = 2) => { if (v == null || !isFinite(v)) return '–'; const r = Number(v.toFixed(d)); return (r < 0 ? MINUS : '+') + Math.abs(r).toFixed(d); };
const pct = (v, d = 0) => (v == null || !isFinite(v)) ? '–' : v.toFixed(d) + '%';

/* ---------- color scale: same diverging ramp and fixed range as the printed maps ---------- */
const STOPS = ['#1b5fa8', '#6fa3d6', '#cfe0f1', '#f0efec', '#f4d3c4', '#e48a6c', '#b8321f'].map((h) => [1, 3, 5].map((i) => parseInt(h.slice(i, i + 2), 16)));
const LUT_N = 512; let VMIN = -2, VMAX = 2; const LUT = new Uint8ClampedArray(LUT_N * 4);
function buildLut() {
  for (let k = 0; k < LUT_N; k++) {
    const t = k / (LUT_N - 1) * (STOPS.length - 1), i = Math.min(Math.floor(t), STOPS.length - 2), f = t - i;
    for (let c = 0; c < 3; c++) LUT[k * 4 + c] = STOPS[i][c] + (STOPS[i + 1][c] - STOPS[i][c]) * f;
    LUT[k * 4 + 3] = 255;
  }
}
const lutIndex = (v) => Math.round((Math.min(Math.max(v, VMIN), VMAX) - VMIN) / (VMAX - VMIN) * (LUT_N - 1));
function colorOf(v) { const k = lutIndex(v) * 4; return [LUT[k], LUT[k + 1], LUT[k + 2]]; }
const inkFor = ([r, g, b]) => (0.299 * r + 0.587 * g + 0.114 * b) < 140 ? '#fff' : '#111';

/* ---------- grid painting ---------- */
function paintTile(canvas, coords, grids) {
  const W = canvas.width, H = canvas.height, S = W / 256, ctx = canvas.getContext('2d'), img = ctx.createImageData(W, H), px = img.data, n = 2 ** coords.z;
  const lons = new Float64Array(W), lats = new Float64Array(H);
  for (let i = 0; i < W; i++) lons[i] = (coords.x + (i + 0.5) / W) / n * 360 - 180;
  for (let j = 0; j < H; j++) lats[j] = Math.atan(Math.sinh(Math.PI * (1 - 2 * (coords.y + (j + 0.5) / H) / n))) * 180 / Math.PI;
  let any = false;
  for (const g of grids) {
    if (lons[0] > g.west + g.nx * g.dx || lons[W - 1] < g.west || lats[0] < g.north - g.ny * g.dy || lats[H - 1] > g.north) continue;
    for (let j = 0; j < H; j++) {
      const iy = Math.floor((g.north - lats[j]) / g.dy); if (iy < 0 || iy >= g.ny) continue;
      const row = iy * g.nx;
      for (let i = 0; i < W; i++) {
        const ix = Math.floor((lons[i] - g.west) / g.dx); if (ix < 0 || ix >= g.nx) continue;
        const v = g.data[row + ix]; if (v !== v) continue;                       // NaN: outside the CWA / too few valid storms
        const o = (j * W + i) * 4; any = true;
        if (v === EXCLUDED) {                                                    // radar beam-blockage sector: gray hatch, anchored to the page grid
          const dark = ((((coords.x * W + i) + (coords.y * H + j)) / S) | 0) % 8 < 2;
          px[o] = dark ? 70 : 170; px[o + 1] = dark ? 70 : 170; px[o + 2] = dark ? 70 : 170; px[o + 3] = dark ? 235 : 190;
        } else { const k = lutIndex(v) * 4; px[o] = LUT[k]; px[o + 1] = LUT[k + 1]; px[o + 2] = LUT[k + 2]; px[o + 3] = 255; }
      }
    }
  }
  if (any) ctx.putImageData(img, 0, 0);
}
const ErrorLayer = L.GridLayer.extend({
  initialize(opts) { L.GridLayer.prototype.initialize.call(this, opts); this.grids = []; },
  createTile(coords) {
    const c = document.createElement('canvas'), S = this.options.scale; c.width = c.height = 256 * S; c.className = 'leaflet-tile';
    paintTile(c, coords, this.grids); return c;
  }
});
function sample(grids, lat, lon) {
  for (const g of grids) {
    const iy = Math.floor((g.north - lat) / g.dy), ix = Math.floor((lon - g.west) / g.dx);
    if (iy < 0 || iy >= g.ny || ix < 0 || ix >= g.nx) continue;
    const v = g.data[iy * g.nx + ix]; if (v === v) return { off: g.off, v };
  }
  return null;
}
function toFloat32LE(buf) {
  if (new Uint8Array(new Uint32Array([1]).buffer)[0] === 1) return new Float32Array(buf);       // little-endian platform: view directly
  const dv = new DataView(buf), out = new Float32Array(buf.byteLength / 4);
  for (let i = 0; i < out.length; i++) out[i] = dv.getFloat32(i * 4, true); return out;
}

/* ---------- state ---------- */
const state = { M: null, season: null, grids: [], sel: null, token: 0, base: 'terrain', opacity: 0.85 };
let map, errLayer, terrainLayer, bm = {}, baseLabelLayer, cwaGeo, cwaCasing, labelLayer, tipBox, legendBox, cwaData;
/* Our own basemap (built by scripts/build_basemap.py): shaded-relief terrain tiles + water / border / state-line vectors + labels, all served from basemap/. */
const BASE_DIR = 'basemap/', LAND_PLAIN = '#eef0ea', WATER_FILL = '#c4d8e4';
const BASEMAPS = { terrain: { name: 'Terrain (shaded relief)' }, plain: { name: 'Plain (no terrain)' } };
const LINE_STYLE = { coast: { color: '#5f7a8a', weight: 1 }, state: { color: '#59626b', weight: 1.1, dashArray: '5 3' }, border: { color: '#2c343b', weight: 2 } };

function notice(msg) { const n = el('notice'); n.hidden = !msg; n.textContent = msg || ''; }
async function getJson(url) { const r = await fetch(url); if (!r.ok) throw new Error(`${url}: HTTP ${r.status}`); return r.json(); }
const seasonData = () => state.M.seasons[state.season];
const visibleGrids = () => state.grids.filter((g) => g.visible);

/* ---------- hash state ---------- */
function readHash() { const p = new URLSearchParams(location.hash.slice(1)); return { season: p.get('season'), office: p.get('office'), base: p.get('base'), op: parseFloat(p.get('op')) }; }
function writeHash() {
  const p = new URLSearchParams({ season: state.season, base: state.base, op: state.opacity.toFixed(2) }); if (state.sel) p.set('office', state.sel);
  history.replaceState(null, '', '#' + p.toString());
}

/* ---------- data loading ---------- */
async function loadSeason(id) {
  const token = ++state.token, s = state.M.seasons[id]; notice('');
  try {
    const grids = await Promise.all(Object.entries(s.offices).map(async ([off, e]) => {
      const r = await fetch(url.bin(id, off, e));
      if (!r.ok) throw new Error(`${off}: HTTP ${r.status}`);
      const buf = await r.arrayBuffer(); if (buf.byteLength !== e.nx * e.ny * 4) throw new Error(`${off}: expected ${e.nx * e.ny * 4} bytes, got ${buf.byteLength}`);
      return Object.assign({ off, visible: true, data: toFloat32LE(buf) }, e);
    }));
    if (token !== state.token) return;                      // a newer season was chosen meanwhile
    el('topSeason').textContent = s.label || id; state.season = id; state.grids = grids; if (state.sel && !s.offices[state.sel]) state.sel = null;
    refresh(); drawLegend();
  } catch (err) { notice('Could not load data: ' + err.message); }
}

/* ---------- rendering ---------- */
function refresh() {
  errLayer.grids = visibleGrids(); errLayer.redraw(); drawTable(); drawLabels(); styleCwa(); drawDetail(); writeHash();
}
function drawTable() {
  const rows = state.grids.map((g) => {
    const st = g.stats, c = colorOf(st.mean_error), cls = ['row', state.sel === g.off ? 'sel' : '', g.visible ? '' : 'off'].join(' '), o = state.M.offices[g.off];
    return `<tr class="${cls}" data-off="${g.off}" tabindex="0">
      <td><input type="checkbox" data-vis="${g.off}" ${g.visible ? 'checked' : ''} aria-label="Show ${esc(g.off)}"></td>
      <td class="codename"><span class="chip" style="background:${esc(o.color)}"></span><b>${esc(g.off)}</b><small>${esc(o.name)}</small></td>
      <td>${st.storms}${st.small_sample ? ' <span class="warn" title="Fewer than ' + state.M.small_sample_storms + ' storms: indicative only">&#9888;</span>' : ''}</td>
      <td><span class="val" style="background:rgb(${c});color:${inkFor(c)}">${fmt(st.mean_error)}</span></td>
      <td>${st.mae.toFixed(2)}</td><td>${pct(st.hit_rate)}</td></tr>`;
  }).join('');
  el('offices').innerHTML = `<thead><tr><th></th><th>Office</th><th>Storms</th><th>Mean error (in)</th><th>MAE (in)</th><th>Hit</th></tr></thead><tbody>${rows}</tbody>`;
}
function drawLabels() {
  labelLayer.clearLayers(); if (!cwaData) return;
  for (const f of cwaData.features) {
    const g = state.grids.find((x) => x.off === f.properties.cwa); if (!g || !g.visible) continue;
    const [lo, la] = f.properties.label;
    labelLayer.addLayer(L.marker([la, lo], { interactive: false, keyboard: false, icon: L.divIcon({ className: '', iconSize: [0, 0], html: `<div class="cwa-label">${esc(g.off)}<span>${fmt(g.stats.mean_error)} in</span></div>` }) }));
  }
}
function styleCwa() {
  if (!cwaGeo) return;
  cwaGeo.eachLayer((l) => { const off = l.feature.properties.cwa, g = state.grids.find((x) => x.off === off), on = !!g && g.visible, sel = state.sel === off;
    l.setStyle({ color: state.M.offices[off].color, weight: sel ? 4.5 : 2.4, opacity: on ? 1 : 0.35, fill: true, fillOpacity: 0.001 }); });
  cwaCasing.eachLayer((l) => { const off = l.feature.properties.cwa, sel = state.sel === off; l.setStyle({ color: '#fff', weight: sel ? 8 : 5, opacity: 0.85, fill: false }); });
}
function drawDetail() {
  const box = el('detail'), g = state.grids.find((x) => x.off === state.sel);
  if (!g) { box.hidden = true; box.innerHTML = ''; return; }
  const st = g.stats, o = state.M.offices[g.off], rx = state.M.radar_exclusions && state.M.radar_exclusions[g.off], ci = st.mean_error_ci90;
  box.hidden = false;
  box.innerHTML = `<h3><span class="chip" style="background:${esc(o.color)}"></span>${esc(g.off)} &mdash; ${esc(o.name)}</h3>
    <div>${esc(seasonData().label)}: ${st.storms} storms, ${esc(st.first_storm)} to ${esc(st.last_storm)}</div>
    ${st.small_sample ? '<p class="warn">Only ' + st.storms + ' storms: treat as indicative, not conclusive.</p>' : ''}
    <dl>
      <dt>Mean error (forecast &minus; observed)</dt><dd>${fmt(st.mean_error)} in</dd>
      <dt>90% storm-bootstrap interval</dt><dd>${fmt(ci[0])} to ${fmt(ci[1])} in</dd>
      <dt>Percent bias</dt><dd>${fmt(st.percent_bias, 1)}%</dd>
      <dt>Mean absolute error</dt><dd>${st.mae.toFixed(2)} in</dd>
      <dt>RMSE</dt><dd>${st.rmse.toFixed(2)} in</dd>
      <dt>Hit rate (within max(2 in, 30% of observed))</dt><dd>${pct(st.hit_rate)}</dd>
      <dt>Mean forecast / mean observed</dt><dd>${st.mean_forecast.toFixed(2)} / ${st.mean_observed.toFixed(2)} in</dd>
      <dt>CWA area with mean error &gt; +0.5 in</dt><dd>${pct(st.area_pct_gt_p05)}</dd>
      <dt>CWA area with mean error &lt; &minus;0.5 in</dt><dd>${pct(st.area_pct_lt_m05)}</dd>
      <dt>Full record (${st.full_record_storms} storms, since 2018)</dt><dd>${fmt(st.full_record_mean_error)} in</dd>
    </dl>
    ${st.cells_excluded_radar ? `<p>${st.cells_excluded_radar} of ${st.cells + st.cells_excluded_radar} CWA cells are hatched and excluded from every statistic: radar beam-blockage artifact${rx ? ' (' + esc(rx.radar) + ', azimuth ' + rx.azimuth[0] + '&deg;&ndash;' + rx.azimuth[1] + '&deg;, range &gt; ' + rx.min_range_km + ' km)' : ''}.</p>` : ''}`;
}
function drawLegend() {
  const anyEx = state.grids.some((g) => g.stats.cells_excluded_radar > 0);
  legendBox.innerHTML = `<div class="ttl">${esc(state.M.layer.label)}</div><div>${esc(state.M.layer.units)}</div><canvas width="270" height="14"></canvas>
    <div class="ticks">${[-2, -1, 0, 1, 2].map((t) => `<span>${t === -2 ? '≤ ' : t === 2 ? '≥ ' : ''}${t === 0 ? '0' : fmt(t, 0)}</span>`).join('')}</div>
    <div class="ends"><span>forecast too low</span><span>forecast too high</span></div>
    ${anyEx ? '<div class="hatch"><i></i>Excluded: radar beam-blockage artifact</div>' : ''}`;
  const c = legendBox.querySelector('canvas'), x = c.getContext('2d');
  for (let i = 0; i < c.width; i++) { const [r, g, b] = colorOf(VMIN + (VMAX - VMIN) * i / (c.width - 1)); x.fillStyle = `rgb(${r},${g},${b})`; x.fillRect(i, 0, 1, c.height); }
}

/* ---------- selection ---------- */
function select(off, zoom = true) {
  state.sel = (state.sel === off && !zoom) ? null : off; styleCwa(); drawTable(); drawDetail(); writeHash();
  const g = state.grids.find((x) => x.off === off);
  if (g && zoom && state.sel) map.fitBounds([[g.north - g.ny * g.dy, g.west], [g.north, g.west + g.nx * g.dx]], { padding: [20, 20], maxZoom: 9 });
}

/* ---------- setup ---------- */
async function loadBasemapFiles() {
  try {
    const [meta, water, lines, labels] = await Promise.all(['meta.json', 'water.geojson', 'lines.geojson', 'labels.json'].map((f) => getJson(BASE_DIR + f)));
    bm = { meta, water, lines, labels };
  } catch (err) { bm = {}; notice('Basemap files missing: run scripts/build_basemap.py (' + err.message + ')'); }
}
function setBasemap(key) {
  if (!BASEMAPS[key]) key = 'terrain'; state.base = key;
  if (terrainLayer) { if (key === 'terrain') { if (!map.hasLayer(terrainLayer)) terrainLayer.addTo(map); } else map.removeLayer(terrainLayer); }
  el('map').style.background = key === 'terrain' ? '#d9dfd0' : LAND_PLAIN; el('basemap').value = key; writeHash();
}
function drawBaseLabels() {
  if (!baseLabelLayer || !bm.labels) return; baseLabelLayer.clearLayers(); const z = map.getZoom();
  for (const l of bm.labels) if (z >= l.minzoom)
    baseLabelLayer.addLayer(L.marker([l.lat, l.lon], { pane: 'baseLabels', interactive: false, keyboard: false, icon: L.divIcon({ className: '', iconSize: [0, 0], html: `<div class="bm-label bm-${l.kind}">${esc(l.name)}</div>` }) }));
}
function initMap() {
  const m = bm.meta, opts = { zoomSnap: 0.25, worldCopyJump: false, attributionControl: false, minZoom: m ? m.minzoom : 5, maxZoom: m ? m.maxzoom + 1 : 10 };
  if (m) { opts.maxBounds = L.latLngBounds(m.bounds); opts.maxBoundsViscosity = 1.0; }
  map = L.map('map', opts).setView([44.6, -71.3], 7);
  if (m) { const fit = () => map.setMinZoom(Math.max(m.minzoom, map.getBoundsZoom(opts.maxBounds, true))); fit(); map.on('resize', fit); }      // never zoom out past the edge of the basemap  // draw order, bottom to top: terrain tiles (tile pane, 200) < water 250 < snowfall-error overlay 300 < borders and coasts 390 < CWA outlines 400 < place names 440 < office labels 600
  for (const [n, z] of [['water', 250], ['data', 300], ['lines', 390], ['baseLabels', 440]]) { const p = map.createPane(n); p.style.zIndex = z; p.style.pointerEvents = 'none'; }
  L.control.attribution({ prefix: false }).addTo(map).addAttribution('Data: NWS NDFD (NCEI), NOHRSC &middot; Terrain: <a href="https://github.com/tilezen/joerd/blob/master/docs/attribution.md" target="_blank" rel="noopener">AWS Terrain Tiles</a> &middot; Boundaries: Natural Earth &middot; CWA: NWS');
  if (m) {
    terrainLayer = L.tileLayer(BASE_DIR + m.tiles, { minNativeZoom: m.minzoom, maxNativeZoom: m.maxzoom, minZoom: m.minzoom, maxZoom: m.maxzoom + 1, bounds: m.bounds, tileSize: 256, zIndex: 1, keepBuffer: 4, errorTileUrl: 'data:image/gif;base64,R0lGODlhAQABAAAAACH5BAEKAAEALAAAAAABAAEAAAICTAEAOw==' });
    L.geoJSON(bm.water, { pane: 'water', renderer: L.svg({ pane: 'water' }), interactive: false, style: { stroke: false, fillColor: WATER_FILL, fillOpacity: 1 } }).addTo(map);
    const lineLayer = (kind, extra) => L.geoJSON(bm.lines, { pane: 'lines', renderer: L.svg({ pane: 'lines' }), interactive: false, filter: (f) => f.properties.kind === kind, style: Object.assign({ lineCap: 'round', lineJoin: 'round' }, LINE_STYLE[kind], extra) }).addTo(map);
    lineLayer('state', { color: '#fff', weight: 3, opacity: 0.55, dashArray: null }); lineLayer('border', { color: '#fff', weight: 4, opacity: 0.55 });      // pale casings keep the lines readable over the data colours
    lineLayer('coast'); lineLayer('state'); lineLayer('border');
    baseLabelLayer = L.layerGroup().addTo(map); drawBaseLabels(); map.on('zoomend', drawBaseLabels);
  }
  errLayer = new ErrorLayer({ pane: 'data', tileSize: 256, scale: Math.max(2, Math.ceil(window.devicePixelRatio || 1)),   /* cells are only ~4 px wide at zoom 7: paint at 2x so fractional zoom levels do not blur them */ zIndex: 5, opacity: state.opacity, updateWhenIdle: false, keepBuffer: 2 }).addTo(map);
  labelLayer = L.layerGroup().addTo(map);
  legendBox = L.control({ position: 'bottomright' }); legendBox.onAdd = () => L.DomUtil.create('div', 'legend'); legendBox.addTo(map); legendBox = legendBox.getContainer();

  // pointer tooltip: follows the cursor over the map
  const tipHtml = (ll) => {
    const s = sample(visibleGrids(), ll.lat, ll.lng), pos = `${ll.lat.toFixed(3)}°N ${Math.abs(ll.lng).toFixed(3)}°W`;
    if (!s) return `<div class="tip-none">No data</div><div class="tip-pos">${pos}</div>`;
    const o = state.M.offices[s.off] || {}, head = `<div class="tip-head"><i class="chip" style="background:${esc(o.color || '#888')}"></i><b>${esc(s.off)}</b>${o.name ? ' <span>' + esc(o.name) + '</span>' : ''}</div>`;
    if (s.v === EXCLUDED) return `${head}<div class="tip-none">Excluded: radar beam-blockage artifact</div><div class="tip-pos">${pos}</div>`;
    const c = colorOf(s.v);
    return `${head}<div class="tip-val"><span class="val" style="background:rgb(${c[0]},${c[1]},${c[2]});color:${inkFor(c)}">${fmt(s.v)} in</span> mean error per storm</div><div class="tip-pos">${pos}</div>`;
  };
  const touch = window.matchMedia('(hover: none)').matches, box = map.getContainer(); let dragging = false, last = null;
  tipBox = L.DomUtil.create('div', 'cursor-tip', box); tipBox.hidden = true;
  const hideTip = () => { tipBox.hidden = true; };
  const showTip = (cp) => {
    last = cp; tipBox.innerHTML = tipHtml(map.containerPointToLatLng(cp)); tipBox.hidden = false;
    const S = map.getSize(), w = tipBox.offsetWidth, h = tipBox.offsetHeight; let x = cp.x + 16, y = cp.y + 18;
    if (x + w > S.x - 6) x = cp.x - w - 12; if (y + h > S.y - 6) y = cp.y - h - 12;
    tipBox.style.transform = `translate(${Math.max(4, Math.round(x))}px, ${Math.max(4, Math.round(y))}px)`;
  };
  if (!touch) {
    box.addEventListener('mousemove', (e) => { if (dragging) return; if (e.target.closest('.leaflet-control, .leaflet-popup')) { last = null; return hideTip(); } showTip(map.mouseEventToContainerPoint(e)); });
    box.addEventListener('mouseleave', () => { last = null; hideTip(); });
    map.on('dragstart', () => { dragging = true; hideTip(); }); map.on('dragend', () => { dragging = false; });
    map.on('zoomstart', hideTip); map.on('moveend', () => { if (last && !dragging) showTip(last); });         // after a zoom or fly-to the cell under a still pointer has changed
    el('offices').addEventListener('change', () => { if (last) showTip(last); });
  } else map.on('click', (e) => L.popup().setLatLng(e.latlng).setContent(tipHtml(e.latlng)).openOn(map));     // touch screens: tap to read a cell
}
async function initCwa() {
  try {
    cwaData = await getJson(url.cwa());
    cwaCasing = L.geoJSON(cwaData, { interactive: false, style: { color: '#fff', weight: 5, fill: false } }).addTo(map);
    cwaGeo = L.geoJSON(cwaData, { style: { weight: 2.4, fill: true, fillOpacity: 0.001 }, onEachFeature: (f, l) => l.on('click', () => select(f.properties.cwa, false)) }).addTo(map);
    cwaGeo.bringToFront(); styleCwa(); drawLabels();
  } catch (err) { notice('Could not load CWA outlines: ' + err.message); }
}

async function main() {
  buildLut();
  const inline = el('manifest'); try { state.M = inline ? JSON.parse(inline.textContent) : await getJson(url.manifest()); } catch (err) { notice('No data found. Run scripts/export_viewer_data.py first.'); return; }
  VMIN = state.M.layer.vmin; VMAX = state.M.layer.vmax; buildLut();
  const h = readHash(), seasons = Object.keys(state.M.seasons).sort().reverse();
  if (isFinite(h.op)) state.opacity = Math.min(1, Math.max(0, h.op));
  el('season').innerHTML = seasons.map((s) => `<option value="${esc(s)}">${esc(state.M.seasons[s].label)}</option>`).join('');
  el('basemap').innerHTML = Object.entries(BASEMAPS).map(([k, b]) => `<option value="${k}">${esc(b.name)}</option>`).join('');
  el('opacity').value = state.opacity; el('opacityOut').textContent = Math.round(state.opacity * 100) + '%';
  await loadBasemapFiles();
  initMap(); window.NWS.map = map; window.NWS.state = state;      // exposed for embedding and testing
  setBasemap(h.base || 'terrain');
  const season = seasons.includes(h.season) ? h.season : (state.M.default_season || seasons[0]); el('season').value = season;
  if (h.office) state.sel = h.office.toUpperCase();
  document.addEventListener('keydown', (e) => { if (e.key === 'Escape' && state.sel) { state.sel = null; styleCwa(); drawTable(); drawDetail(); writeHash(); } });
  el('season').addEventListener('change', (e) => loadSeason(e.target.value));
  el('basemap').addEventListener('change', (e) => setBasemap(e.target.value));
  el('opacity').addEventListener('input', (e) => { state.opacity = parseFloat(e.target.value); errLayer.setOpacity(state.opacity); el('opacityOut').textContent = Math.round(state.opacity * 100) + '%'; writeHash(); });
  el('offices').addEventListener('change', (e) => { const off = e.target.dataset && e.target.dataset.vis; if (!off) return; state.grids.find((g) => g.off === off).visible = e.target.checked; refresh(); });
  el('offices').addEventListener('click', (e) => { if (e.target.closest('input')) return; const tr = e.target.closest('tr.row'); if (tr) select(tr.dataset.off, true); });
  el('offices').addEventListener('keydown', (e) => { if (e.key === 'Enter' && e.target.matches('tr.row')) select(e.target.dataset.off, true); });
  await Promise.all([loadSeason(season), initCwa()]);
  if (state.grids.length) {
    if (state.sel && state.grids.some((g) => g.off === state.sel)) select(state.sel, true);
    else map.fitBounds([[Math.min(...state.grids.map((g) => g.north - g.ny * g.dy)), Math.min(...state.grids.map((g) => g.west))],
                        [Math.max(...state.grids.map((g) => g.north)), Math.max(...state.grids.map((g) => g.west + g.nx * g.dx))]], { padding: [10, 10] });
  }
  drawLabels(); styleCwa();
}
main();
})();
