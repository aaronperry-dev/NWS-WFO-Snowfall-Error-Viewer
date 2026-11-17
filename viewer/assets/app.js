/* NWS WFO snowfall forecast error viewer.
 * Loads one raw float32 .bin grid per office (via api.php), and paints them on a Leaflet map with a custom tile layer that samples the
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
const state = { M: null, season: null, grids: [], sel: null, token: 0, base: 'light', opacity: 0.85 };
let map, errLayer, baseLayer, labelsTiles, cwaGeo, cwaCasing, labelLayer, readoutBox, legendBox, cwaData;
const BASEMAPS = {
  light: { name: 'Light (CARTO)', url: 'https://{s}.basemaps.cartocdn.com/light_all/{z}/{x}/{y}{r}.png', o: { subdomains: 'abcd', maxZoom: 19, attribution: '&copy; OpenStreetMap contributors &copy; CARTO' } },
  streets: { name: 'Streets (OpenStreetMap)', url: 'https://tile.openstreetmap.org/{z}/{x}/{y}.png', o: { maxZoom: 19, attribution: '&copy; OpenStreetMap contributors' } },
  topo: { name: 'Topographic (Esri)', url: 'https://server.arcgisonline.com/ArcGIS/rest/services/World_Topo_Map/MapServer/tile/{z}/{y}/{x}', o: { maxZoom: 19, attribution: 'Tiles &copy; Esri' } },
  imagery: { name: 'Satellite (Esri)', url: 'https://server.arcgisonline.com/ArcGIS/rest/services/World_Imagery/MapServer/tile/{z}/{y}/{x}', o: { maxZoom: 19, attribution: 'Tiles &copy; Esri' } }
};

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
    state.season = id; state.grids = grids; if (state.sel && !s.offices[state.sel]) state.sel = null;
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
function setBasemap(key) {
  if (!BASEMAPS[key]) key = 'light'; state.base = key; if (baseLayer) map.removeLayer(baseLayer);
  baseLayer = L.tileLayer(BASEMAPS[key].url, Object.assign({ zIndex: 1 }, BASEMAPS[key].o)).addTo(map); el('basemap').value = key; writeHash();
  if (labelsTiles) { map.removeLayer(labelsTiles); labelsTiles = null; }
  if (key === 'light') labelsTiles = L.tileLayer('https://{s}.basemaps.cartocdn.com/light_only_labels/{z}/{x}/{y}{r}.png', { pane: 'basemapLabels', subdomains: 'abcd', maxZoom: 19 }).addTo(map);   // place names stay readable above the data
}
function initMap() {
  map = L.map('map', { zoomSnap: 0.25, worldCopyJump: false, attributionControl: false }).setView([44.6, -71.3], 7);
  map.createPane('basemapLabels').style.zIndex = 450; map.getPane('basemapLabels').style.pointerEvents = 'none';
  L.control.attribution({ prefix: false }).addTo(map).addAttribution('Data: NWS NDFD (NCEI), NOHRSC &middot; CWA: NWS');
  errLayer = new ErrorLayer({ tileSize: 256, scale: Math.max(2, Math.ceil(window.devicePixelRatio || 1)),   /* cells are only ~4 px wide at zoom 7: paint at 2x so fractional zoom levels do not blur them */ zIndex: 5, opacity: state.opacity, updateWhenIdle: false, keepBuffer: 2 }).addTo(map);
  labelLayer = L.layerGroup().addTo(map);
  readoutBox = L.control({ position: 'bottomleft' }); readoutBox.onAdd = () => { const d = L.DomUtil.create('div', 'readout'); d.innerHTML = 'Move the pointer over the map'; return d; }; readoutBox.addTo(map);
  legendBox = L.control({ position: 'bottomright' }); legendBox.onAdd = () => L.DomUtil.create('div', 'legend'); legendBox.addTo(map); legendBox = legendBox.getContainer();
  const readoutHtml = (ll) => {
    const s = sample(visibleGrids(), ll.lat, ll.lng), pos = `${ll.lat.toFixed(3)}°N ${Math.abs(ll.lng).toFixed(3)}°W`;
    return !s ? `<span style="color:#5b6670">No data here</span><br>${pos}`
      : s.v === EXCLUDED ? `<b>${esc(s.off)}</b> &middot; excluded (radar artifact)<br>${pos}`
      : `<b>${esc(s.off)}</b> &middot; mean error <b>${fmt(s.v)} in</b><br>${pos}`;
  };
  map.on('mousemove', (e) => { readoutBox.getContainer().innerHTML = readoutHtml(e.latlng); });
  if (window.matchMedia('(hover: none)').matches) map.on('click', (e) => L.popup().setLatLng(e.latlng).setContent(readoutHtml(e.latlng)).openOn(map));     // touch screens: tap to read a cell
  map.on('mouseout', () => { readoutBox.getContainer().innerHTML = 'Move the pointer over the map'; });
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
  initMap(); window.NWS.map = map; window.NWS.state = state;      // exposed for embedding and testing
  setBasemap(h.base || 'light');
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
