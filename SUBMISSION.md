# Space Apps Global Submission — ready to paste

**Challenge:** SPHEREx (Challenge 12)
**Project name:** SPHEREx Agent Hunt
**Link to Final Project:** https://github.com/mrfentmen/spherex-agent-hunt
**Demo:** 7 slides — `slides.html` in the repo (open in any browser)
**License:** Apache 2.0
**AI use disclosure:** Super Intelligence

## Summary (short)

A crew of AI agents hunts for undiscovered moving objects in NASA's SPHEREx infrared sky survey. Nothing gets called a discovery until it survives six proof gates and an independent Skeptic re-run. The human is not inside the hunt; the human is the last pair of eyes who receives the news only when a candidate survives everything.

## Description

SPHEREx has photographed the whole sky in 102 infrared bands. Asteroids, comets, and distant slow movers shift a little between visits; stars sit still and glitches lie. No human can eyeball patch after patch, epoch after epoch, so the hunt is fully agentic.

Seven agents, one job each: Scout picks sky patches, Fetcher pulls SPHEREx images from IRSA over the IVOA SIA v2 protocol (no key, no login), Watcher groups images into epochs, Tracker links detections and fits straight-line motion, Classifier labels tracks, Skeptic tries to kill every candidate on sabotaged control data, and Reporter writes the case file with full provenance (IRSA file IDs and access links for every claim).

Six proof gates decide: at least 3 repeat epochs, linear motion, a hold-out prediction test, artifact exclusion, catalog cross-match against Gaia / WISE / MPC, and spectrum consistency across the 102 bands. A FAIL kills a candidate. UNVERIFIED keeps it a candidate. A case file requires all six to PASS.

Calibration, on synthetic data with a planted mover (labelled test data): the hunt finds the mover at 32.33 arcsec/day across 4 epochs, fit rms 0.09 arcsec, hold-out prediction off by 0.26 arcsec against a 6.0 tolerance. Stars are killed by the motion gate, glitches by roundness, and no case file is issued without a catalog. 10 of 10 unit tests pass.

Live data status, stated plainly: live IRSA inventory is verified working (10 real spherex_qr3 images pulled through the pipeline's own client on Nov 15, 2026, after a client bug was found by bisection and fixed). No real-sky hunt has run yet, and no discovery is claimed. Not yet built: FITS source extraction (detection tables are the current input) and automated catalog fetch (catalogs are CSV-supplied for now).

## Data and sources

- SPHEREx image metadata via IRSA, the NASA/IPAC Infrared Science Archive (SIA v2 endpoint, collections spherex_dr1 and spherex_qr3)
- IVOA Simple Image Access protocol, version 2
- Gaia, WISE, and Minor Planet Center catalogs (cross-match targets)
- Python standard library only; no third-party dependencies
