"""Fetcher layer: query IRSA for SPHEREx images via IVOA SIA v2.

Verified against IRSA's SPHEREx data access documentation
(https://caltech-ipac.github.io/spherex-archive-documentation/spherex-data-access/):

- IRSA's generic SIA2 endpoint is https://irsa.ipac.caltech.edu/SIA
- SPHEREx collections searchable via SIA2 include spherex_dr1 and
  spherex_qr3 (DR1 and QR3 supersede QR2; QR2 is not searchable via SIA2).
- A query returns a VOTable whose rows carry an access_url per image.
- No API key, account, or login is required.

Only the Python standard library is used, so this runs anywhere.
"""
from __future__ import annotations

import csv
import json
import urllib.parse
import urllib.request
import xml.etree.ElementTree as ET
from dataclasses import asdict, dataclass

SIA_ENDPOINT = "https://irsa.ipac.caltech.edu/SIA"
COLLECTIONS = ("spherex_dr1", "spherex_qr3")
USER_AGENT = "spherex-agent-hunt/0.1 (NASA Space Apps 2026)"


@dataclass
class ImageRecord:
    """One SPHEREx image from an SIA query, with provenance."""

    obs_id: str
    collection: str
    ra_deg: float | None
    dec_deg: float | None
    t_min_mjd: float | None
    t_max_mjd: float | None
    access_url: str
    instrument: str = ""
    dataproduct_type: str = ""

    def to_dict(self) -> dict:
        return asdict(self)


def build_query_url(ra_deg: float, dec_deg: float, radius_deg: float,
                    collection: str, maxrec: int = 2000) -> str:
    if collection not in COLLECTIONS:
        raise ValueError(f"unknown SPHEREx collection {collection!r}; "
                         f"expected one of {COLLECTIONS}")
    params = {
        "POS": f"CIRCLE {ra_deg} {dec_deg} {radius_deg}",
        "COLLECTION": collection,
        "MAXREC": str(maxrec),
        # No FORMAT parameter: IRSA's SIA service silently returns an
        # empty table (QUERY_STATUS OK, zero rows) when
        # FORMAT=application/x-votable+xml is passed. Verified
        # 2026-11-15 by bisection: same query without FORMAT returns
        # the real rows; the default response is already VOTable.
    }
    return SIA_ENDPOINT + "?" + urllib.parse.urlencode(params)


def _local(tag: str) -> str:
    return tag.rsplit("}", 1)[-1]


def parse_votable(xml_bytes: bytes) -> list[ImageRecord]:
    """Parse an ObsCore VOTable returned by IRSA's SIA service."""
    root = ET.fromstring(xml_bytes)
    table = None
    for el in root.iter():
        if _local(el.tag) == "TABLE":
            table = el
            break
    if table is None:
        raise ValueError("no TABLE element in VOTable response")
    fields: list[str] = []
    for child in table:
        if _local(child.tag) == "FIELD":
            fields.append(child.get("name", ""))
    records: list[ImageRecord] = []
    for tr in table.iter():
        if _local(tr.tag) != "TR":
            continue
        cells = []
        for td in tr:
            if _local(td.tag) == "TD":
                cells.append((td.text or "").strip())
        row = dict(zip(fields, cells))

        def fnum(*names: str) -> float | None:
            for n in names:
                v = row.get(n, "")
                if v:
                    try:
                        return float(v)
                    except ValueError:
                        return None
            return None

        records.append(ImageRecord(
            obs_id=row.get("obs_id", "") or row.get("obs_publisher_did", ""),
            collection=row.get("obs_collection", ""),
            ra_deg=fnum("s_ra", "ra"),
            dec_deg=fnum("s_dec", "dec"),
            t_min_mjd=fnum("t_min"),
            t_max_mjd=fnum("t_max"),
            access_url=row.get("access_url", ""),
            instrument=row.get("instrument_name", ""),
            dataproduct_type=row.get("dataproduct_type", ""),
        ))
    return records


def query_sia(ra_deg: float, dec_deg: float, radius_deg: float,
              collection: str = "spherex_qr3", maxrec: int = 2000,
              timeout: float = 300.0) -> list[ImageRecord]:
    """Run a live SIA2 cone query against IRSA. Returns image records."""
    url = build_query_url(ra_deg, dec_deg, radius_deg, collection, maxrec)
    req = urllib.request.Request(url, headers={"User-Agent": USER_AGENT})
    with urllib.request.urlopen(req, timeout=timeout) as resp:
        payload = resp.read()
    return parse_votable(payload)


def write_inventory(records: list[ImageRecord], path: str) -> None:
    if path.endswith(".json"):
        with open(path, "w") as fh:
            json.dump([r.to_dict() for r in records], fh, indent=2)
        return
    cols = ["obs_id", "collection", "ra_deg", "dec_deg", "t_min_mjd",
            "t_max_mjd", "access_url", "instrument", "dataproduct_type"]
    with open(path, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=cols)
        w.writeheader()
        for r in records:
            w.writerow(r.to_dict())
