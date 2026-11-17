<?php
declare(strict_types=1);
/**
 * Viewer configuration.
 * DATA_DIR holds manifest.json, cwa.geojson and <season>/<OFFICE>_mean_error.bin, as written by scripts/export_viewer_data.py.
 * Override with the NWS_VIEWER_DATA environment variable to keep the data outside the web root.
 */
const APP_TITLE = 'NWS WFO Snowfall Forecast Error Viewer';
$dataDir = getenv('NWS_VIEWER_DATA') ?: __DIR__ . '/data';
define('DATA_DIR', rtrim($dataDir, '/'));

/** Parsed manifest.json, or null if the exporter has not been run yet. */
function load_manifest(): ?array {
    $f = DATA_DIR . '/manifest.json';
    if (!is_file($f)) return null;
    $m = json_decode((string) file_get_contents($f), true);
    return is_array($m) ? $m : null;
}
