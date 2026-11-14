# SPHEREx Agent Hunt — Sunday pitch notes

**When:** Sun Nov 15, 2026 — slides due 12:00 PM, presentations start 12:30 PM, **5 min + 2 min Q&A**, 370 Jay St, Brooklyn.
**Deck:** `slides.html` in this repo — open in a browser, present full-screen (or print to PDF from the browser). 7 slides, ~40 seconds each.
**Rule:** everything on the slides is verified from the repo. The calibration numbers are synthetic test data and the slides say so. No discovery is claimed anywhere.

## Talk track (5 minutes)

**1. Title (~20s)** — "SPHEREx Agent Hunt. A crew of AI agents hunting for undiscovered moving objects in NASA's SPHEREx infrared survey — and nothing gets called a discovery until it survives every proof check."

**2. Problem (~40s)** — SPHEREx photographed the whole sky in 102 infrared bands. Asteroids, comets and distant slow movers shift a little between visits; stars sit still, glitches lie. One human can't eyeball patch after patch, epoch after epoch — so the hunt is fully agentic and the human is only the last pair of eyes.

**3. Crew (~45s)** — Seven agents, one job each: Scout picks patches, Fetcher pulls SPHEREx images from IRSA, Watcher groups epochs, Tracker fits straight-line motion, Classifier labels tracks, Skeptic tries to kill every candidate on sabotaged data, Reporter writes the case file with provenance.

**4. Proof bar (~45s)** — Six gates: repeat epochs (3+), linear motion, hold-out prediction, artifact exclusion, catalog cross-match (Gaia/WISE/MPC), spectrum consistency across 102 bands. FAIL kills a candidate. UNVERIFIED keeps it a candidate. A case file needs all six PASS.

**5. Calibration receipt (~60s)** — This is synthetic test data, labelled on the slide. Planted mover: found at 32.33 arcsec/day, 4 epochs, fit rms 0.09 arcsec, hold-out prediction off by 0.26 arcsec. Six stars killed as stationary, two glitches killed on roundness. Zero case files — the mover stays a candidate until the catalog gate passes. 10/10 tests pass.

**6. Live data (~40s)** — Honest status: the IRSA client works and parses real VOTable replies, no login needed. But live inventory from the build sandbox failed all day Saturday (HTTP 502 proxy error); it answered fine in prep with 428 images. Next: retry from a normal connection, run a real patch hunt, add FITS extraction and automated catalog fetch.

**7. Close (~20s)** — Repo is public, Apache 2.0, stdlib-only Python. AI use disclosed as "Super Intelligence". "Try it. Check it. Kill it if it's wrong."

## Likely Q&A

- **"Did you find anything real?"** — No. No real-sky hunt has run yet; the calibration mover is synthetic. The pipeline is proven, the sky run is next.
- **"Why should we trust an agent's discovery?"** — You shouldn't trust the agent. That's the point of the six gates plus the Skeptic: the system has to survive its own sabotage controls before a human ever sees it.
- **"What breaks first on real data?"** — Catalog cross-match and spectrum consistency are UNVERIFIED until real catalogs and per-band fluxes are supplied. FITS source extraction is not built yet.
- **"Solo team?"** — Yes, one person. The agents are the crew; the human is the last pair of eyes.

## Numbers source (re-verified Sat Nov 14, 13:43 ET)

- Tests: `PYTHONPATH=src python3 -m unittest discover -s tests` — 10/10 pass.
- Calibration hunt: `calibrate` → 44 rows, 4 epochs; `hunt` → planted mover 32.33 arcsec/day, rms 0.09, hold-out 0.26 arcsec (tolerance 6.0); stars 0.00–0.01 arcsec/day killed; glitches killed on roundness 0.9; 0 case files.
- Live IRSA retry same afternoon: HTTP 502 Proxy Error in ~1.6s — still blocked from this sandbox, not a code failure (client parsed a valid empty VOTable earlier).
- Skeptic caveat, in the report: time-shuffle control dissolved tracks (8 vs 9); position-scramble still linked 9 tracks — flagged as WARNING in the report, open work, not hidden.
