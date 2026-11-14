"""Classifier + Skeptic layer: the proof gates.

Del's rule for this project: nothing gets called a discovery until it
survives every proof check. Each gate returns PASS, FAIL, or UNVERIFIED:

- FAIL kills the candidate outright.
- UNVERIFIED means the data needed for that check does not exist yet.
  The candidate stays a candidate. It is never promoted on a promise.
- A case file (the only output a human is asked to look at) requires
  every required gate to PASS.

Gates, from the project spec:
1. repeat_epochs     - detected in 3 or more distinct epochs
2. linear_motion     - one straight-line motion fits all epochs,
                       residuals small, rate in a physical range
3. holdout_prediction- fit all epochs but the last; the fit must predict
                       the last epoch's position within tolerance
4. artifact_exclusion- shape metrics look like a point source, not a
                       glitch: roundness, pixel area, SNR per epoch
5. catalog_crossmatch- no known stationary catalog source at the track's
                       mean position (Gaia / WISE / MPC exports, CSV).
                       UNVERIFIED when no catalog file was supplied.
6. spectrum_consistency - SPHEREx sees 102 spectral bands; a real source
                       keeps a consistent per-band flux pattern across
                       epochs. UNVERIFIED when band fluxes are absent.
"""
from __future__ import annotations

import csv
import math
from dataclasses import dataclass, field

from .detections import Track, angular_sep_arcsec, fit_motion

PASS, FAIL, UNVERIFIED = "PASS", "FAIL", "UNVERIFIED"

# SPHEREx pixel scale is about 6.2 arcsec/pixel; a point source spans a
# few pixels. Tolerances below are in those physical terms.
DEFAULTS = {
    "min_epochs": 3,
    "max_rms_arcsec": 3.0,
    "min_rate_arcsec_day": 1.0,      # slower than this is a star, not a mover
    "max_rate_arcsec_day": 1800.0,   # 0.5 deg/day ceiling for a real track
    "holdout_tolerance_arcsec": 6.0,
    "max_roundness": 0.5,
    "min_pixel_area": 2.0,
    "max_pixel_area": 200.0,
    "min_snr": 5.0,
    "catalog_radius_arcsec": 3.0,
    "max_band_scatter": 0.35,        # fractional scatter across epochs
}


@dataclass
class GateResult:
    gate: str
    status: str
    detail: str


@dataclass
class Verdict:
    gates: list[GateResult] = field(default_factory=list)

    @property
    def is_case_file(self) -> bool:
        return bool(self.gates) and all(g.status == PASS for g in self.gates)

    @property
    def is_killed(self) -> bool:
        return any(g.status == FAIL for g in self.gates)

    @property
    def label(self) -> str:
        if self.is_case_file:
            return "CASE FILE - survives every proof gate"
        if self.is_killed:
            return "KILLED - failed a proof gate"
        return "CANDIDATE - one or more gates still unverified"


def _gate_epochs(track: Track, cfg: dict) -> GateResult:
    n = track.epochs
    ok = n >= cfg["min_epochs"]
    return GateResult("repeat_epochs", PASS if ok else FAIL,
                      f"{n} distinct epochs (need {cfg['min_epochs']})")


def _gate_linear(track: Track, cfg: dict) -> GateResult:
    fit = track.fit or fit_motion(track.detections)
    if fit.rate_arcsec_per_day < cfg["min_rate_arcsec_day"]:
        return GateResult("linear_motion", FAIL,
                          f"rate {fit.rate_arcsec_per_day:.2f} arcsec/day is "
                          f"stationary, not a mover")
    if fit.rate_arcsec_per_day > cfg["max_rate_arcsec_day"]:
        return GateResult("linear_motion", FAIL,
                          f"rate {fit.rate_arcsec_per_day:.1f} arcsec/day is "
                          f"above the physical ceiling")
    if fit.rms_arcsec > cfg["max_rms_arcsec"]:
        return GateResult("linear_motion", FAIL,
                          f"fit rms {fit.rms_arcsec:.2f} arcsec exceeds "
                          f"{cfg['max_rms_arcsec']} arcsec")
    return GateResult("linear_motion", PASS,
                      f"rate {fit.rate_arcsec_per_day:.2f} arcsec/day, "
                      f"rms {fit.rms_arcsec:.2f} arcsec over "
                      f"{track.span_days:.1f} days")


def _gate_holdout(track: Track, cfg: dict) -> GateResult:
    dets = sorted(track.detections, key=lambda d: d.mjd)
    if len({d.epoch_id for d in dets}) < 3 or len(dets) < 3:
        return GateResult("holdout_prediction", UNVERIFIED,
                          "need at least 3 detections to hold one out")
    train, held = dets[:-1], dets[-1]
    try:
        fit = fit_motion(train)
    except ValueError as exc:
        return GateResult("holdout_prediction", UNVERIFIED, str(exc))
    pra, pdec = fit.predict(held.mjd)
    dist = angular_sep_arcsec(pra, pdec, held.ra_deg, held.dec_deg)
    ok = dist <= cfg["holdout_tolerance_arcsec"]
    return GateResult("holdout_prediction", PASS if ok else FAIL,
                      f"predicted the held-out epoch within {dist:.2f} "
                      f"arcsec (tolerance {cfg['holdout_tolerance_arcsec']})")


def _gate_artifact(track: Track, cfg: dict) -> GateResult:
    for d in track.detections:
        if d.roundness > cfg["max_roundness"]:
            return GateResult("artifact_exclusion", FAIL,
                              f"epoch {d.epoch_id}: roundness {d.roundness} "
                              f"is elongated, glitch-like")
        if not (cfg["min_pixel_area"] <= d.pixel_area <= cfg["max_pixel_area"]):
            return GateResult("artifact_exclusion", FAIL,
                              f"epoch {d.epoch_id}: pixel area {d.pixel_area} "
                              f"outside point-source range")
        if d.snr < cfg["min_snr"]:
            return GateResult("artifact_exclusion", FAIL,
                              f"epoch {d.epoch_id}: SNR {d.snr} below "
                              f"{cfg['min_snr']}")
    return GateResult("artifact_exclusion", PASS,
                      "every epoch is round, point-sized, and above SNR floor")


def load_catalog(path: str) -> list[tuple[float, float, str]]:
    """Catalog rows: ra_deg, dec_deg, name (Gaia/WISE/MPC export, CSV)."""
    rows = []
    with open(path, newline="") as fh:
        for row in csv.DictReader(fh):
            rows.append((float(row["ra_deg"]), float(row["dec_deg"]),
                         row.get("name", "catalog source")))
    return rows


def _gate_catalog(track: Track, cfg: dict,
                  catalog: list[tuple[float, float, str]] | None) -> GateResult:
    if catalog is None:
        return GateResult("catalog_crossmatch", UNVERIFIED,
                          "no known-source catalog supplied yet "
                          "(Gaia / WISE / MPC)")
    fit = track.fit or fit_motion(track.detections)
    for ra, dec, name in catalog:
        dist = angular_sep_arcsec(fit.ra0, fit.dec0, ra, dec)
        if dist <= cfg["catalog_radius_arcsec"]:
            return GateResult("catalog_crossmatch", FAIL,
                              f"matches known source {name} at {dist:.2f} arcsec")
    return GateResult("catalog_crossmatch", PASS,
                      f"no known source within {cfg['catalog_radius_arcsec']} "
                      f"arcsec of the track position")


def _gate_spectrum(track: Track, cfg: dict) -> GateResult:
    by_band: dict[str, list[float]] = {}
    for d in track.detections:
        if d.band and d.flux_band is not None:
            by_band.setdefault(d.band, []).append(d.flux_band)
    usable = {b: v for b, v in by_band.items() if len(v) >= 2}
    if not usable:
        return GateResult("spectrum_consistency", UNVERIFIED,
                          "no per-band fluxes across epochs yet "
                          "(SPHEREx 102-band check pending)")
    for band, vals in usable.items():
        mean = sum(vals) / len(vals)
        if mean <= 0:
            return GateResult("spectrum_consistency", FAIL,
                              f"band {band}: non-positive mean flux")
        scatter = math.sqrt(sum((v - mean) ** 2 for v in vals) / len(vals)) / mean
        if scatter > cfg["max_band_scatter"]:
            return GateResult("spectrum_consistency", FAIL,
                              f"band {band}: flux scatter {scatter:.2f} across "
                              f"epochs is not a consistent spectrum")
    return GateResult("spectrum_consistency", PASS,
                      f"per-band flux consistent across epochs in "
                      f"{len(usable)} band(s)")


def evaluate(track: Track, catalog=None, cfg: dict | None = None) -> Verdict:
    cfg = {**DEFAULTS, **(cfg or {})}
    return Verdict(gates=[
        _gate_epochs(track, cfg),
        _gate_linear(track, cfg),
        _gate_holdout(track, cfg),
        _gate_artifact(track, cfg),
        _gate_catalog(track, cfg, catalog),
        _gate_spectrum(track, cfg),
    ])


def classify(track: Track, verdict: Verdict) -> str:
    """Plain-language class for one track."""
    fit = track.fit
    if fit and fit.rate_arcsec_per_day < DEFAULTS["min_rate_arcsec_day"]:
        return "stationary source (star-like)"
    if verdict.is_killed:
        failed = [g.gate for g in verdict.gates if g.status == FAIL]
        return f"rejected ({', '.join(failed)})"
    if verdict.is_case_file:
        return "moving-object candidate, all proof gates passed"
    return "moving-object candidate, proof incomplete"
