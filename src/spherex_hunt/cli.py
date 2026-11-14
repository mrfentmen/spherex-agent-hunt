"""Command line entry point.

    python -m spherex_hunt inventory --ra 210.80225 --dec 54.34894 \
        --radius 0.1 --collection spherex_qr3 --out inventory.csv

    python -m spherex_hunt calibrate --out calibration.csv

    python -m spherex_hunt hunt --name patch-1 --ra 210.80225 \
        --dec 54.34894 --detections calibration.csv \
        [--catalog known_sources.csv] --out runs/patch-1
"""
from __future__ import annotations

import argparse
import os
import sys

from . import agents, calibration, sia
from .proofs import load_catalog


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(prog="spherex_hunt")
    sub = ap.add_subparsers(dest="cmd", required=True)

    inv = sub.add_parser("inventory", help="live IRSA SIA query for a patch")
    inv.add_argument("--ra", type=float, required=True)
    inv.add_argument("--dec", type=float, required=True)
    inv.add_argument("--radius", type=float, default=0.1)
    inv.add_argument("--collection", default="spherex_qr3",
                     choices=sia.COLLECTIONS)
    inv.add_argument("--out", required=True)

    cal = sub.add_parser("calibrate", help="write synthetic calibration data")
    cal.add_argument("--out", required=True)

    hunt = sub.add_parser("hunt", help="run the crew over a detection table")
    hunt.add_argument("--name", required=True)
    hunt.add_argument("--ra", type=float, required=True)
    hunt.add_argument("--dec", type=float, required=True)
    hunt.add_argument("--radius", type=float, default=0.25)
    hunt.add_argument("--detections", required=True)
    hunt.add_argument("--catalog")
    hunt.add_argument("--out", required=True)

    args = ap.parse_args(argv)

    if args.cmd == "inventory":
        records = sia.query_sia(args.ra, args.dec, args.radius,
                                collection=args.collection)
        sia.write_inventory(records, args.out)
        print(f"{len(records)} images -> {args.out}")
        epochs = agents.watcher(records)
        print(f"epochs by time: { {k: len(v) for k, v in epochs.items()} }")
        return 0

    if args.cmd == "calibrate":
        info = calibration.generate(args.out)
        print(f"calibration data (synthetic): {info} -> {args.out}")
        return 0

    if args.cmd == "hunt":
        patch = agents.Patch(args.name, args.ra, args.dec, args.radius)
        if args.detections.endswith(".csv") and _is_calibration(args.detections):
            detections = calibration.load_calibration(args.detections)
        else:
            from .detections import load_detections
            detections = load_detections(args.detections)
        catalog = load_catalog(args.catalog) if args.catalog else None
        result = agents.run_hunt(patch, detections, catalog=catalog)
        os.makedirs(args.out, exist_ok=True)
        report = agents.reporter(result)
        with open(os.path.join(args.out, "report.txt"), "w") as fh:
            fh.write(report + "\n")
        n_case = 0
        for i, (_, verdict, _) in enumerate(result.tracks):
            if verdict.is_case_file:
                n_case += 1
                with open(os.path.join(args.out, f"case-file-{i + 1}.md"),
                          "w") as fh:
                    fh.write(agents.case_file_markdown(result, i) + "\n")
        print(report)
        print(f"report + {n_case} case file(s) -> {args.out}")
        return 0
    return 2


def _is_calibration(path: str) -> bool:
    with open(path) as fh:
        return fh.readline().startswith("# CALIBRATION DATA")


if __name__ == "__main__":
    sys.exit(main())
