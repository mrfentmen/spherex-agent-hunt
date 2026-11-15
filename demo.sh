#!/bin/bash
# Offline demo: proves the pipeline with zero network. Exit nonzero on any failure.
set -euo pipefail
cd "$(dirname "$0")"
echo "== tests =="
PYTHONPATH=src python3 -m unittest discover -s tests 2>&1 | tail -3
echo "== calibration hunt (synthetic) =="
PYTHONPATH=src python3 -m spherex_hunt calibrate --out /tmp/demo_cal.csv
PYTHONPATH=src python3 -m spherex_hunt hunt --name demo --ra 210.80225 --dec 54.34894 \
  --detections /tmp/demo_cal.csv --out /tmp/runs/demo 2>&1 | grep -E "Tracks linked|Case files"
echo "demo ok: no network touched, nothing claimed"
