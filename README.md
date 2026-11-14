# SPHEREx Agent Hunt

**NASA Space Apps Challenge 2026, Challenge 12: Planet X and SPHEREx.**
Team: Abdelrahman El Sayad (solo), New York local event.

A fully agentic hunter for undiscovered moving objects (asteroids,
comets, and distant slow movers) in NASA's SPHEREx infrared sky survey.
A crew of agents hunts on its own. The human is only the last pair of
eyes, and only sees a candidate after it survives every proof gate.

## The rule

Nothing gets called a discovery until it survives every proof check.
Each gate returns PASS, FAIL, or UNVERIFIED. FAIL kills a candidate.
UNVERIFIED keeps it a candidate, never a discovery. A **case file**,
the only output a human is asked to review, requires all gates to PASS.

## The crew

| Agent | Job |
|---|---|
| Scout | Picks the sky patches to hunt, in priority order |
| Fetcher | Queries IRSA for SPHEREx images covering a patch (IVOA SIA v2) |
| Watcher | Groups a patch's images into observation epochs by time |
| Tracker | Links per-epoch detections into straight-line motion tracks |
| Classifier | Labels tracks: stationary, moving candidate, rejected |
| Skeptic | Re-runs the hunt on sabotaged data (shuffled times, scrambled positions). Real tracks must dissolve under sabotage |
| Reporter | Writes the run report and any case files, with full provenance |

## Proof gates

1. **Repeat epochs**: detected in 3 or more distinct epochs.
2. **Linear motion**: one straight-line fit covers all epochs, small
   residuals, rate in a physical range (stationary = star, killed).
3. **Hold-out prediction**: fit all epochs but the last; the fit must
   predict the held-out epoch's position within tolerance.
4. **Artifact exclusion**: every epoch is round, point-sized, and above
   the SNR floor (kills cosmic rays and glitches).
5. **Catalog cross-match**: no known Gaia / WISE / MPC source at the
   track position. UNVERIFIED until a catalog export is supplied.
6. **Spectrum consistency**: SPHEREx observes in 102 spectral bands; a
   real source keeps a consistent per-band flux across epochs.
   UNVERIFIED until per-band fluxes are supplied.

## Quick start

No dependencies. Python 3.10+, standard library only.

```bash
# Live: ask IRSA what SPHEREx images cover a patch of sky
PYTHONPATH=src python -m spherex_hunt inventory \
  --ra 210.80225 --dec 54.34894 --radius 0.1 \
  --collection spherex_qr3 --out inventory.csv

# Calibration: synthetic data with a planted mover (TEST DATA ONLY)
PYTHONPATH=src python -m spherex_hunt calibrate --out calibration.csv

# Hunt: link detections, run the proof gates, write the report
PYTHONPATH=src python -m spherex_hunt hunt \
  --name patch-1 --ra 210.80225 --dec 54.34894 \
  --detections calibration.csv --out runs/patch-1

# Tests (ground truth: the planted mover must be found, stars and
# glitches must not become case files)
PYTHONPATH=src python -m unittest discover -s tests -v
```

## Data and provenance

- SPHEREx images are served by **IRSA** (NASA/IPAC Infrared Science
  Archive) through the IVOA Simple Image Access v2 service at
  `https://irsa.ipac.caltech.edu/SIA`, collections `spherex_dr1` and
  `spherex_qr3`. No API key, account, or login is required. Every image
  record keeps its IRSA `obs_id` and `access_url`, and every detection
  row carries the `image_id` it came from, so each claim in a case file
  traces back to a specific IRSA file.
- Detection tables are CSV with columns: `source_id, epoch_id, mjd,
  ra_deg, dec_deg, snr, flux, roundness, pixel_area, x_px, y_px,
  image_id, band, flux_band`.

## Status, stated plainly

- **Working and tested**: SIA v2 query client and VOTable parsing;
  epoch grouping; track linking and motion fitting; all six proof
  gates; Skeptic negative controls; case-file and report generation;
  calibration data with known ground truth.
- **Not yet built**: automatic source extraction from downloaded SPHEREx
  FITS images (needs astropy/SEP; detection tables are currently the
  pipeline's input), and an automated Gaia/WISE/MPC catalog fetch
  (catalogs are currently supplied as CSV exports). These are labeled
  UNVERIFIED in the gates rather than silently assumed.

## AI use disclosure

This project was built with the assistance of **Super Intelligence**
(artificial intelligence), per the Space Apps AI-use disclosure rules.
All data are NASA data; all findings require human verification.

## References

- NASA Space Apps Challenge 2026, Challenge 12: Planet X and SPHEREx.
  https://www.spaceappschallenge.org
- SPHEREx Data Access, IRSA:
  https://caltech-ipac.github.io/spherex-archive-documentation/spherex-data-access/
- IRSA SIA v2 service documentation:
  https://irsa.ipac.caltech.edu/ibe/sia.html
- IRSA SPHEREx cutout tutorial (worked SIA/TAP examples):
  https://caltech-ipac.github.io/irsa-tutorials/spherex-cutouts/
- IVOA Simple Image Access v2 protocol: https://www.ivoa.net/documents/SIA/

## License

Apache License 2.0. See [LICENSE](LICENSE).
