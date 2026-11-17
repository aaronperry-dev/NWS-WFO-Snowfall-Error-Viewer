<?php
declare(strict_types=1);
require __DIR__ . '/config.php';
$manifest = load_manifest();
$v = fn(string $f): string => (string) @filemtime(__DIR__ . '/' . $f);      // cache-busting version for assets
$h = fn(string $s): string => htmlspecialchars($s, ENT_QUOTES | ENT_SUBSTITUTE, 'UTF-8');
?><!doctype html>
<html lang="en">
<head>
<meta charset="utf-8">
<meta name="viewport" content="width=device-width, initial-scale=1">
<title><?= $h(APP_TITLE) ?></title>
<link rel="stylesheet" href="vendor/leaflet/leaflet.css">
<link rel="stylesheet" href="assets/style.css?v=<?= $v('assets/style.css') ?>">
</head>
<body>
<header class="topbar">
  <h1><?= $h(APP_TITLE) ?></h1>
  <p>New England NWS forecast offices &middot; <span id="topSeason"></span></p>
</header>
<div class="app">
  <div id="map" role="application" aria-label="Map of mean snowfall forecast error by NWS office"></div>
  <aside class="panel">
    <header>
      <p class="sub">Mean snowfall forecast error (forecast − observed) for the New England NWS forecast offices, drawn at the native ~4&nbsp;km grid of the observations.</p>
    </header>
    <div id="notice" class="notice" hidden></div>
    <section class="controls">
      <label>Season <select id="season"></select></label>
      <label>Basemap <select id="basemap"></select></label>
      <label>Overlay opacity <input id="opacity" type="range" min="0" max="1" step="0.05" value="0.85"> <output id="opacityOut">85%</output></label>
    </section>
    <section>
      <h2>Offices</h2>
      <table id="offices" class="offices"></table>
      <p class="hint">Tick to show or hide an office. Click a row (or a CWA on the map) for details.</p>
    </section>
    <section id="detail" class="detail" hidden></section>
    <section class="about">
      <h2>How to read this</h2>
      <ul>
        <li><b>Red</b> = the forecast was too high on average; <b>blue</b> = too low. Each cell is the mean of (forecast − observed) storm-total snowfall over the season's storms, in inches per storm.</li>
        <li>Forecast: the official WFO-edited NDFD snowfall grid, issued about 12&nbsp;h before each storm. Observed: NOHRSC National Snowfall Analysis. Storms: CWA-mean total of 1.5&nbsp;in or more.</li>
        <li>Single-season averages rest on few storms (5&ndash;24 per office); intervals are wide. See each office's storm count and 90&nbsp;% interval.</li>
        <li>Hatched cells are excluded radar beam-blockage artifacts in the observation analysis.</li>
      </ul>
    </section>
  </aside>
</div>
<?php if ($manifest !== null): ?>
<script id="manifest" type="application/json"><?= json_encode($manifest, JSON_HEX_TAG | JSON_HEX_AMP | JSON_HEX_APOS | JSON_HEX_QUOT | JSON_UNESCAPED_SLASHES) ?></script>
<?php endif; ?>
<script>window.NWS = { api: 'api.php' };</script>
<script src="vendor/leaflet/leaflet.js"></script>
<script src="assets/app.js?v=<?= $v('assets/app.js') ?>"></script>
<noscript><p style="padding:1rem">This viewer needs JavaScript.</p></noscript>
</body>
</html>
