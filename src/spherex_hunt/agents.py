"""The agent crew, as plain functions with one job each.

Scout    - picks the sky patches to hunt (from a target list).
Fetcher  - asks IRSA what SPHEREx images cover a patch (SIA v2).
Watcher  - groups a patch's images into observation epochs.
Tracker  - links per-epoch detections into straight-line tracks.
Classifier - labels each track (see proofs.classify).
Skeptic  - re-runs the hunt on sabotaged data: shuffled epoch times and
           an offset empty-sky position. A real pipeline must find the
           candidate in the real data and NOT find it in the sabotage.
Reporter - writes the case file for any track that survives every gate,
           and a run report that lists kills and open candidates too.

The human is the last pair of eyes: only case files are presented as
discoveries, and even they are labeled as candidates for human review.
"""
from __future__ import annotations

import random
from dataclasses import dataclass, field

from . import sia
from .detections import (Detection, Track, angular_sep_arcsec, fit_motion,
                         link_tracks)
from .proofs import Verdict, classify, evaluate


@dataclass
class Patch:
    name: str
    ra_deg: float
    dec_deg: float
    radius_deg: float = 0.25


def scout(patches: list[Patch]) -> list[Patch]:
    """Scout: today the target list is the strategy. Order is priority."""
    return list(patches)


def fetcher(patch: Patch, collection: str = "spherex_qr3") -> list[sia.ImageRecord]:
    """Fetcher: live IRSA SIA v2 query for one patch."""
    return sia.query_sia(patch.ra_deg, patch.dec_deg, patch.radius_deg,
                         collection=collection)


def watcher(records: list[sia.ImageRecord],
            gap_days: float = 0.5) -> dict[str, list[sia.ImageRecord]]:
    """Watcher: group images into epochs by observation time."""
    timed = sorted((r for r in records if r.t_min_mjd is not None),
                   key=lambda r: r.t_min_mjd)
    epochs: dict[str, list[sia.ImageRecord]] = {}
    current: list[sia.ImageRecord] = []
    last_t = None
    for r in timed:
        if last_t is not None and r.t_min_mjd - last_t > gap_days:
            epochs[f"epoch_{len(epochs) + 1}"] = current
            current = []
        current.append(r)
        last_t = r.t_min_mjd
    if current:
        epochs[f"epoch_{len(epochs) + 1}"] = current
    return epochs


def tracker(detections: list[Detection], **kw) -> list[Track]:
    return link_tracks(detections, **kw)


@dataclass
class HuntResult:
    patch: Patch
    images: list[sia.ImageRecord] = field(default_factory=list)
    tracks: list[tuple[Track, Verdict, str]] = field(default_factory=list)
    skeptic_notes: list[str] = field(default_factory=list)


def skeptic(detections: list[Detection], real_track_count: int,
            seed: int = 3301) -> list[str]:
    """Skeptic negative controls. Returns human-readable notes.

    Control 1, time shuffle: permute the MJDs among detections. Straight
    lines in the real data must dissolve; if shuffled data links just as
    many tracks, the linker is manufacturing candidates.
    Control 2, position scramble: rotate all positions by a random sky
    offset per detection. Coherent tracks must dissolve.
    """
    notes: list[str] = []
    rng = random.Random(seed)
    if len(detections) >= 3:
        times = [d.mjd for d in detections]
        shuffled_times = times[:]
        rng.shuffle(shuffled_times)
        sabotaged = [Detection(**{**d.__dict__, "mjd": t})
                     for d, t in zip(detections, shuffled_times)]
        n = len(link_tracks(sabotaged))
        notes.append(
            f"time-shuffle control: {n} track(s) from shuffled times vs "
            f"{real_track_count} real. "
            + ("OK, sabotage dissolved the tracks." if n < max(real_track_count, 1)
               else "WARNING: shuffled data links as well as real data."))
        scrambled = [Detection(**{**d.__dict__,
                                  "ra_deg": d.ra_deg + rng.uniform(-0.01, 0.01),
                                  "dec_deg": d.dec_deg + rng.uniform(-0.01, 0.01)})
                     for d in detections]
        n2 = len(link_tracks(scrambled))
        notes.append(
            f"position-scramble control: {n2} track(s) from scrambled "
            f"positions. "
            + ("OK." if n2 < max(real_track_count, 1) else
               "WARNING: scrambled data still links."))
    else:
        notes.append("skeptic controls skipped: fewer than 3 detections")
    return notes


def run_hunt(patch: Patch, detections: list[Detection],
             catalog=None, fetch_live: bool = False) -> HuntResult:
    """Run the full crew over one patch's detection table."""
    result = HuntResult(patch=patch)
    if fetch_live:
        result.images = fetcher(patch)
    tracks = tracker(detections)
    for tr in tracks:
        if tr.fit is None:
            tr.fit = fit_motion(tr.detections)
        verdict = evaluate(tr, catalog=catalog)
        result.tracks.append((tr, verdict, classify(tr, verdict)))
    result.skeptic_notes = skeptic(detections, len(tracks))
    return result


def reporter(result: HuntResult) -> str:
    """Reporter: the run report, in plain text a tired human can scan."""
    p = result.patch
    lines = [
        f"SPHEREx hunt report: {p.name} "
        f"(RA {p.ra_deg}, Dec {p.dec_deg}, radius {p.radius_deg} deg)",
        f"Images in inventory: {len(result.images)}",
        f"Tracks linked: {len(result.tracks)}",
        "",
    ]
    for i, (track, verdict, label) in enumerate(result.tracks, 1):
        fit = track.fit
        lines.append(f"Track {i}: {label}")
        lines.append(f"  verdict: {verdict.label}")
        if fit:
            lines.append(
                f"  motion: {fit.rate_arcsec_per_day:.2f} arcsec/day, "
                f"rms {fit.rms_arcsec:.2f} arcsec, "
                f"{track.epochs} epochs over {track.span_days:.1f} days")
        for g in verdict.gates:
            lines.append(f"  gate {g.gate}: {g.status} - {g.detail}")
        lines.append("  provenance (IRSA image ids):")
        seen = []
        for d in track.detections:
            if d.image_id and d.image_id not in seen:
                seen.append(d.image_id)
        lines.extend(f"    {img}" for img in seen or ["(none recorded)"])
        lines.append("")
    lines.append("Skeptic negative controls:")
    lines.extend(f"  {n}" for n in result.skeptic_notes)
    case_files = [v for _, v, _ in result.tracks if v.is_case_file]
    lines.append("")
    lines.append(f"Case files for human review: {len(case_files)}")
    lines.append("Nothing in this report is a confirmed discovery. "
                 "Case files are candidates that survived every gate.")
    return "\n".join(lines)


def case_file_markdown(result: HuntResult, index: int) -> str:
    """One candidate's case file: the only document a human must read."""
    track, verdict, label = result.tracks[index]
    fit = track.fit
    lines = [
        f"# Case file: {result.patch.name} track {index + 1}",
        "",
        f"Class: {label}",
        f"Verdict: {verdict.label}",
        "",
        "## Motion",
        f"- Rate: {fit.rate_arcsec_per_day:.2f} arcsec/day",
        f"- Fit RMS: {fit.rms_arcsec:.2f} arcsec",
        f"- Epochs: {track.epochs} over {track.span_days:.1f} days",
        f"- Position at fit epoch (MJD {fit.t_ref:.1f}): "
        f"RA {fit.ra0:.5f}, Dec {fit.dec0:.5f}",
        "",
        "## Proof gates",
    ]
    for g in verdict.gates:
        lines.append(f"- {g.gate}: {g.status} - {g.detail}")
    lines += ["", "## Detections (provenance)"]
    for d in sorted(track.detections, key=lambda x: x.mjd):
        lines.append(
            f"- MJD {d.mjd:.2f} RA {d.ra_deg:.5f} Dec {d.dec_deg:.5f} "
            f"SNR {d.snr} roundness {d.roundness} area {d.pixel_area} px "
            f"image {d.image_id or '(unrecorded)'}")
    lines += ["", "## Skeptic controls", ""]
    lines.extend(f"- {n}" for n in result.skeptic_notes)
    return "\n".join(lines)
