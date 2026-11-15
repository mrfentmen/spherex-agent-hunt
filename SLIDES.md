# SPHEREx Agent Hunt — Sunday pitch notes

**When:** Sun Nov 15, 2026 — slides due 12:00 PM, presentations start 12:30 PM, **5 min + 2 min Q&A**, 370 Jay St, Brooklyn.
**Deck:** `slides.html` in this repo — open in a browser, present full-screen (or print to PDF from the browser). 7 slides, ~40 seconds each.
**Rule:** everything on the slides is verified from the repo. The calibration numbers are synthetic test data and the slides say so. No discovery is claimed anywhere.

## Talk track (5 minutes)

**1. Title (~20s)** — "SPHEREx Agent Hunt. A crew of AI agents hunting for undiscovered moving objects in NASA's SPHEREx infrared survey — and nothing gets called a discovery until it survives every proof check."

**2. Problem (~40s)** — SPHEREx photographed the whole sky in 102 infrared bands. Asteroids, comets and distant slow movers shift a little between visits; stars sit still, glitches lie. One human can't eyeball patch after patch, epoch after epoch — so the hunt is fully agentic and the human is only the last pair of eyes.

**3. Crew (~45s)** — Seven agents, one job each: Scout picks patches, Fetcher pulls SPHEREx images from IRSA, Watcher groups epochs, Tracker fits straight-line motion, Classifier labels tracks, Skeptic tries to kill every candidate on sabotaged data, Reporter writes the case file with provenance.

**4. Proof bar (~45s)** — Six gates: repeat epochs (3+), linear motion, hold-out prediction, artifact exclusion, catalog cross-match (Gaia/WISE/MPC), spectrum consistency across 102 bands. FAIL kills a candidate. UNVERIFIED keeps it a candidate. A case file needs all six PASS.

**5. Calibration receipt (~60s)** — This is synthetic test data, labelled on the slide. Planted mover: found at 32.33 arcsec/day, 4 epochs, fit rms 0.09 arcsec, hold-out prediction off by 0.26 arcsec. Six stars killed as stationary, two glitches killed on roundness. Zero case files — the mover stays a candidate until the catalog gate passes. 64/64 tests pass.

**6. Live data (~40s)** — Live inventory VERIFIED Sun Nov 15: the empty results were our own bug, a FORMAT parameter IRSA silently answers with an empty table; found by bisection, fixed, pushed — the same query then returned real rows. A 0.05 deg cone returns 132 SPHEREx images in 15 epochs, no login. Then a real 3-frame hunt across 36 days: 180 detections, 2 tracks linked through all 3 epochs — and the artifact gate killed both. Zero case files. Wide cones still die at the sandbox proxy, so hunts run in small patches.

**7. Close (~20s)** — Repo is public, Apache 2.0, stdlib-only Python. AI use disclosed as "Super Intelligence". "Try it. Check it. Kill it if it's wrong."

## Likely Q&A

- **"Did you find anything real?"** — Two real tracks linked across 3 epochs and 36 days — and the artifact gate killed both. Zero case files. The calibration mover is synthetic; the pipeline plus the kills are proven on both.
- **"Why should we trust an agent's discovery?"** — You shouldn't trust the agent. That's the point of the six gates plus the Skeptic: the system has to survive its own sabotage controls before a human ever sees it.
- **"What breaks first on real data?"** — Artifact shapes and the time-shuffle control: both real tracks died on shape, and shuffled times still linked 2 tracks, flagged WARNING. Catalog cross-match stays UNVERIFIED until real catalogs are supplied.
- **"Solo team?"** — Yes, one person. The agents are the crew; the human is the last pair of eyes.

## Numbers source (re-verified Sun Nov 15, morning ET)

- Tests: `PYTHONPATH=src python3 -m unittest discover -s tests` — 64/64 pass.
- Calibration hunt: `calibrate` → 44 rows, 4 epochs; `hunt` → planted mover 32.33 arcsec/day, rms 0.09, hold-out 0.26 arcsec (tolerance 6.0); stars 0.00–0.01 arcsec/day killed; glitches killed on roundness 0.9; 0 case files.
- Live IRSA: root cause of the empty results found and fixed by bisection — the client sent FORMAT=application/x-votable+xml, which IRSA silently answers with an empty table; dropping it made the same query return real rows. 0.05 deg cone on spherex_qr3 → 132 images, 15 epochs (same cone on dr1 → 254 images, 26 epochs). Wide queries still die at the sandbox proxy (HTTP 502); small patches pass.
- Real 3-frame hunt: 3 overlapping frames, 36-day span, 180 numpy detections → 2 tracks linked across all 3 epochs (114 and 234 arcsec/day, hold-out PASS), both KILLED on artifact_exclusion, 0 case files. Time-shuffle control: 2 vs 2, flagged NOTE. Full-SIP astrometry check: same-detector pair disagrees by ~500 arcsec bulk, ~17 arcsec residual after alignment — the residual is the accuracy floor.
