<?php
declare(strict_types=1);
/**
 * Read-only data API for the viewer.
 *   api.php?action=manifest                         manifest.json
 *   api.php?action=cwa                              simplified CWA outlines (GeoJSON)
 *   api.php?action=bin&season=2025-26&office=CAR    raw float32 grid for one office and season
 * The season and office are looked up as KEYS of the manifest, never used to build a path, so only files the exporter listed can be served.
 * Responses carry an ETag (304 when unchanged), are gzip-compressed when the client allows it, and are cacheable for an hour.
 */
require __DIR__ . '/config.php';

function fail(int $code, string $msg): never {
    http_response_code($code); header('Content-Type: application/json; charset=utf-8'); header('X-Content-Type-Options: nosniff');
    echo json_encode(['error' => $msg]); exit;
}

function send_file(string $path, string $type): never {
    $real = realpath($path); $root = realpath(DATA_DIR);
    if ($real === false || $root === false || !str_starts_with($real, $root . DIRECTORY_SEPARATOR) || !is_file($real)) fail(404, 'file not found');
    $etag = '"' . md5($real . filemtime($real) . filesize($real)) . '"';
    header('ETag: ' . $etag); header('Cache-Control: public, max-age=3600'); header('X-Content-Type-Options: nosniff'); header('Vary: Accept-Encoding');
    if (($_SERVER['HTTP_IF_NONE_MATCH'] ?? '') === $etag) { http_response_code(304); exit; }
    header('Content-Type: ' . $type);
    $body = (string) file_get_contents($real);
    if (function_exists('gzencode') && str_contains($_SERVER['HTTP_ACCEPT_ENCODING'] ?? '', 'gzip')) { $body = gzencode($body, 6); header('Content-Encoding: gzip'); }
    header('Content-Length: ' . strlen($body));
    if ($_SERVER['REQUEST_METHOD'] !== 'HEAD') echo $body;
    exit;
}

if (!in_array($_SERVER['REQUEST_METHOD'], ['GET', 'HEAD'], true)) { header('Allow: GET, HEAD'); fail(405, 'method not allowed'); }
$m = load_manifest();
if ($m === null) fail(503, 'no data: run scripts/export_viewer_data.py');

switch ($_GET['action'] ?? '') {
    case 'manifest':
        send_file(DATA_DIR . '/manifest.json', 'application/json; charset=utf-8');
    case 'cwa':
        send_file(DATA_DIR . '/cwa.geojson', 'application/geo+json; charset=utf-8');
    case 'bin':
        $season = (string) ($_GET['season'] ?? ''); $office = strtoupper((string) ($_GET['office'] ?? ''));
        $entry = $m['seasons'][$season]['offices'][$office] ?? null;
        if (!is_array($entry) || !isset($entry['file'])) fail(404, 'unknown season or office');
        send_file(DATA_DIR . '/' . $entry['file'], 'application/octet-stream');
    default:
        fail(400, 'unknown action');
}
