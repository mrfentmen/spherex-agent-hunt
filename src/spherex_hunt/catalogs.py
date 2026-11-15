"""Keyless catalog fetch: Gaia DR3 cone via TAP, written as the hunt's catalog CSV.

No dependency beyond the standard library. TAP responses are plain CSV,
so parsing is testable offline; only fetch_catalog needs the network.
"""
from __future__ import annotations

import csv
import io
import urllib.parse
import urllib.request

GAIA_TAP_SYNC = "https://gea.esac.esa.int/tap-server/tap/sync"
USER_AGENT = "spherex-agent-hunt/0.1 (NASA Space Apps 2026)"


def parse_tap_csv(text: str) -> list[tuple[float, float, str]]:
    """TAP CSV (ra, dec, source_id columns) -> [(ra, dec, name)]."""
    rows: list[tuple[float, float, str]] = []
    reader = csv.DictReader(io.StringIO(text))
    for row in reader:
        try:
            ra, dec = float(row["ra"]), float(row["dec"])
        except (KeyError, TypeError, ValueError):
            continue
        sid = (row.get("source_id") or "").strip()
        rows.append((ra, dec, f"Gaia DR3 {sid}" if sid else "Gaia DR3 source"))
    return rows


def fetch_gaia_cone(ra_deg: float, dec_deg: float, radius_deg: float,
                    limit: int = 2000, timeout: float = 120.0) -> list[tuple[float, float, str]]:
    """Live Gaia DR3 cone search. Returns [(ra, dec, name)]."""
    adql = ("SELECT TOP %d ra, dec, phot_g_mean_mag, source_id FROM gaiadr3.gaia_source "
            "WHERE 1=CONTAINS(POINT('ICRS',ra,dec), CIRCLE('ICRS',%f,%f,%f))" % (
                limit, ra_deg, dec_deg, radius_deg))
    params = urllib.parse.urlencode(
        {"REQUEST": "doQuery", "LANG": "ADQL", "FORMAT": "csv", "QUERY": adql})
    req = urllib.request.Request(GAIA_TAP_SYNC + "?" + params,
                                 headers={"User-Agent": USER_AGENT})
    with urllib.request.urlopen(req, timeout=timeout) as resp:
        body = resp.read().decode("utf-8", "replace")
    return parse_tap_csv(body)


def write_catalog_csv(rows: list[tuple[float, float, str]], path: str) -> None:
    with open(path, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=["ra_deg", "dec_deg", "name"])
        w.writeheader()
        for ra, dec, name in rows:
            w.writerow({"ra_deg": ra, "dec_deg": dec, "name": name})
