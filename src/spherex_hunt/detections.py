"""Watcher + Tracker layer: link per-epoch source detections into tracks.

Input is a detection table (CSV) with one row per source detection:

    source_id, epoch_id, mjd, ra_deg, dec_deg, snr, flux,
    roundness, pixel_area, x_px, y_px, image_id

image_id is the IRSA obs_id of the image the detection came from, so every
detection keeps its provenance. Optional columns: band, flux_band.

Linking rule (the same idea real asteroid linking pipelines use):
detections in different epochs belong to one track when a single
straight-line motion in RA/Dec vs time fits all of them within a
positional tolerance. Stationary sources are the zero-motion special
case and are classified as such, never as movers.
"""
from __future__ import annotations

import csv
import math
from dataclasses import dataclass, field

ARCSEC_PER_DEG = 3600.0


@dataclass
class Detection:
    source_id: str
    epoch_id: str
    mjd: float
    ra_deg: float
    dec_deg: float
    snr: float
    flux: float
    roundness: float = 0.0      # 0 = perfectly round, 1 = very elongated
    pixel_area: float = 0.0
    x_px: float = 0.0
    y_px: float = 0.0
    image_id: str = ""
    band: str = ""
    flux_band: float | None = None


def load_detections(path: str) -> list[Detection]:
    out: list[Detection] = []
    with open(path, newline="") as fh:
        for row in csv.DictReader(fh):
            def f(name: str, default: float = 0.0) -> float:
                v = (row.get(name) or "").strip()
                return float(v) if v else default

            fb = (row.get("flux_band") or "").strip()
            out.append(Detection(
                source_id=row.get("source_id", ""),
                epoch_id=row.get("epoch_id", ""),
                mjd=f("mjd"),
                ra_deg=f("ra_deg"),
                dec_deg=f("dec_deg"),
                snr=f("snr"),
                flux=f("flux"),
                roundness=f("roundness"),
                pixel_area=f("pixel_area"),
                x_px=f("x_px"),
                y_px=f("y_px"),
                image_id=row.get("image_id", ""),
                band=row.get("band", ""),
                flux_band=float(fb) if fb else None,
            ))
    return out


def angular_sep_arcsec(a_ra: float, a_dec: float, b_ra: float, b_dec: float) -> float:
    dra = (a_ra - b_ra) * math.cos(math.radians((a_dec + b_dec) / 2.0))
    ddec = a_dec - b_dec
    return math.hypot(dra, ddec) * ARCSEC_PER_DEG


@dataclass
class MotionFit:
    """Least-squares straight-line fit of position vs time."""

    ra0: float          # RA at t_ref (deg)
    dec0: float         # Dec at t_ref (deg)
    vra: float          # deg/day in RA (already includes cos(dec) factor space)
    vdec: float         # deg/day in Dec
    t_ref: float        # MJD
    rms_arcsec: float
    max_resid_arcsec: float
    rate_arcsec_per_day: float

    def predict(self, mjd: float) -> tuple[float, float]:
        dt = mjd - self.t_ref
        return self.ra0 + self.vra * dt, self.dec0 + self.vdec * dt


def fit_motion(dets: list[Detection]) -> MotionFit:
    if len(dets) < 2:
        raise ValueError("need at least 2 detections to fit motion")
    t_ref = sum(d.mjd for d in dets) / len(dets)
    cosd = math.cos(math.radians(sum(d.dec_deg for d in dets) / len(dets)))

    def linfit(ts: list[float], xs: list[float]) -> tuple[float, float]:
        n = len(ts)
        mt = sum(ts) / n
        mx = sum(xs) / n
        sxx = sum((t - mt) ** 2 for t in ts)
        if sxx == 0:
            return mx, 0.0
        sxy = sum((t - mt) * (x - mx) for t, x in zip(ts, xs))
        slope = sxy / sxx
        return mx - slope * mt, slope

    ts = [d.mjd - t_ref for d in dets]
    # Work in a local tangent plane: x = dRA * cos(dec), y = dDec, both in deg.
    ra_ref = dets[0].ra_deg
    dec_ref = dets[0].dec_deg
    xs = [(d.ra_deg - ra_ref) * cosd for d in dets]
    ys = [d.dec_deg - dec_ref for d in dets]
    x0, vx = linfit(ts, xs)
    y0, vy = linfit(ts, ys)
    residuals = []
    for d, t in zip(dets, ts):
        px, py = x0 + vx * t, y0 + vy * t
        ox, oy = (d.ra_deg - ra_ref) * cosd, d.dec_deg - dec_ref
        residuals.append(math.hypot(px - ox, py - oy) * ARCSEC_PER_DEG)
    rms = math.sqrt(sum(r * r for r in residuals) / len(residuals))
    # Convert tangent-plane intercept/slope back to sky coordinates.
    return MotionFit(
        ra0=ra_ref + x0 / cosd,
        dec0=dec_ref + y0,
        vra=vx / cosd,
        vdec=vy,
        t_ref=t_ref,
        rms_arcsec=rms,
        max_resid_arcsec=max(residuals),
        rate_arcsec_per_day=math.hypot(vx, vy) * ARCSEC_PER_DEG,
    )


@dataclass
class Track:
    detections: list[Detection] = field(default_factory=list)
    fit: MotionFit | None = None

    @property
    def epochs(self) -> int:
        return len({d.epoch_id for d in self.detections})

    @property
    def span_days(self) -> float:
        ts = [d.mjd for d in self.detections]
        return max(ts) - min(ts) if len(ts) > 1 else 0.0


def link_tracks(detections: list[Detection],
                tolerance_arcsec: float = 6.0,
                min_epochs: int = 3) -> list[Track]:
    """Greedy tracklet linking.

    Seeds a track from detections in two different epochs, extends it with
    any detection in a new epoch that the pair's motion predicts within
    tolerance, refits, and keeps tracks reaching min_epochs. Each
    detection is used by at most one returned track.
    """
    by_epoch: dict[str, list[Detection]] = {}
    for d in detections:
        by_epoch.setdefault(d.epoch_id, []).append(d)
    epoch_ids = sorted(by_epoch, key=lambda e: min(d.mjd for d in by_epoch[e]))
    if len(epoch_ids) < 2:
        return []
    used: set[int] = set()
    tracks: list[Track] = []
    first, rest = epoch_ids[0], epoch_ids[1:]
    for seed in by_epoch[first]:
        for partner_ep in rest:
            for partner in by_epoch[partner_ep]:
                if id(seed) in used or id(partner) in used:
                    continue
                pair = [seed, partner]
                try:
                    fit = fit_motion(pair)
                except ValueError:
                    continue
                # A pair with an absurd implied rate is not a seed.
                if fit.rate_arcsec_per_day > 3600.0:
                    continue
                members = list(pair)
                member_epochs = {seed.epoch_id, partner.epoch_id}
                for ep in epoch_ids:
                    if ep in member_epochs:
                        continue
                    best, best_dist = None, tolerance_arcsec
                    for cand in by_epoch[ep]:
                        if id(cand) in used or cand in members:
                            continue
                        pra, pdec = fit.predict(cand.mjd)
                        dist = angular_sep_arcsec(pra, pdec, cand.ra_deg, cand.dec_deg)
                        if dist < best_dist:
                            best, best_dist = cand, dist
                    if best is not None:
                        members.append(best)
                        member_epochs.add(ep)
                        fit = fit_motion(members)
                if len(member_epochs) >= min_epochs:
                    members.sort(key=lambda d: d.mjd)
                    track = Track(detections=members, fit=fit_motion(members))
                    tracks.append(track)
                    used.update(id(d) for d in members)
                    break
            if id(seed) in used:
                break
    return tracks
