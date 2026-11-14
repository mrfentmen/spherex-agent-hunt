"""End-to-end tests on calibration data with known ground truth.

Run: python -m unittest discover -s tests -v   (from the repo root, with
PYTHONPATH=src) or pip install -e . first.
"""
import os
import sys
import tempfile
import unittest

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "src"))

from spherex_hunt import agents, calibration, proofs, sia  # noqa: E402
from spherex_hunt.detections import Detection, fit_motion, link_tracks  # noqa: E402

VOTABLE = b"""<?xml version="1.0"?>
<VOTABLE xmlns="http://www.ivoa.net/xml/VOTable/v1.3">
 <RESOURCE><TABLE>
  <FIELD name="obs_id"/><FIELD name="obs_collection"/>
  <FIELD name="s_ra"/><FIELD name="s_dec"/>
  <FIELD name="t_min"/><FIELD name="t_max"/>
  <FIELD name="access_url"/><FIELD name="instrument_name"/>
  <FIELD name="dataproduct_type"/>
  <DATA><TABLEDATA>
   <TR><TD>img-1</TD><TD>spherex_qr3</TD><TD>210.8</TD><TD>54.35</TD>
       <TD>60800.0</TD><TD>60800.01</TD>
       <TD>https://example.invalid/img-1.fits</TD><TD>SPHEREx</TD>
       <TD>image</TD></TR>
   <TR><TD>img-2</TD><TD>spherex_qr3</TD><TD>210.81</TD><TD>54.36</TD>
       <TD>60802.0</TD><TD>60802.01</TD>
       <TD>https://example.invalid/img-2.fits</TD><TD>SPHEREx</TD>
       <TD>image</TD></TR>
  </TABLEDATA></DATA>
 </TABLE></RESOURCE>
</VOTABLE>"""


class TestSiaParsing(unittest.TestCase):
    def test_parse_votable(self):
        records = sia.parse_votable(VOTABLE)
        self.assertEqual(len(records), 2)
        self.assertEqual(records[0].obs_id, "img-1")
        self.assertEqual(records[0].collection, "spherex_qr3")
        self.assertAlmostEqual(records[1].t_min_mjd, 60802.0)
        self.assertTrue(records[0].access_url.endswith("img-1.fits"))

    def test_query_url_rejects_unknown_collection(self):
        with self.assertRaises(ValueError):
            sia.build_query_url(1.0, 2.0, 0.1, "spherex_nope")


class TestCalibrationHunt(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.tmp = tempfile.mkdtemp()
        cls.csv_path = os.path.join(cls.tmp, "calibration.csv")
        calibration.generate(cls.csv_path)
        cls.detections = calibration.load_calibration(cls.csv_path)
        patch = agents.Patch("calibration-patch", 210.80225, 54.34894, 0.25)
        cls.result = agents.run_hunt(patch, cls.detections)

    def test_planted_mover_is_found(self):
        mover_tracks = [
            t for t, _, _ in self.result.tracks
            if any(d.source_id == "PLANTED-MOVER" for d in t.detections)
        ]
        self.assertEqual(len(mover_tracks), 1)
        fit = mover_tracks[0].fit
        # Planted rate was sqrt(30^2 + 12^2) ~ 32.3 arcsec/day.
        self.assertAlmostEqual(fit.rate_arcsec_per_day, 32.3, delta=2.0)
        self.assertGreaterEqual(mover_tracks[0].epochs, 4)

    def test_stars_and_glitches_do_not_become_case_files(self):
        for track, verdict, _ in self.result.tracks:
            ids = {d.source_id for d in track.detections}
            if ids & {"PLANTED-MOVER"}:
                continue
            self.assertFalse(verdict.is_case_file,
                             f"unexpected case file from {ids}")

    def test_mover_gates_without_catalog_stay_unverified(self):
        for track, verdict, _ in self.result.tracks:
            if any(d.source_id == "PLANTED-MOVER" for d in track.detections):
                by_gate = {g.gate: g.status for g in verdict.gates}
                self.assertEqual(by_gate["repeat_epochs"], proofs.PASS)
                self.assertEqual(by_gate["linear_motion"], proofs.PASS)
                self.assertEqual(by_gate["holdout_prediction"], proofs.PASS)
                self.assertEqual(by_gate["artifact_exclusion"], proofs.PASS)
                # Honest labeling: no catalog supplied, so this gate must
                # stay UNVERIFIED and the track must NOT be a case file.
                self.assertEqual(by_gate["catalog_crossmatch"],
                                 proofs.UNVERIFIED)
                self.assertFalse(verdict.is_case_file)

    def test_mover_becomes_case_file_with_clean_catalog(self):
        catalog = [(100.0, -20.0, "far-away source")]
        for track, _, _ in self.result.tracks:
            if any(d.source_id == "PLANTED-MOVER" for d in track.detections):
                verdict = proofs.evaluate(track, catalog=catalog)
                by_gate = {g.gate: g.status for g in verdict.gates}
                self.assertEqual(by_gate["catalog_crossmatch"], proofs.PASS)
                self.assertEqual(by_gate["spectrum_consistency"], proofs.PASS)
                self.assertTrue(verdict.is_case_file)

    def test_catalog_match_kills(self):
        track = next(t for t, _, _ in self.result.tracks
                     if any(d.source_id == "PLANTED-MOVER"
                            for d in t.detections))
        fit = track.fit
        verdict = proofs.evaluate(
            track, catalog=[(fit.ra0, fit.dec0, "KNOWN-STAR-1")])
        by_gate = {g.gate: g.status for g in verdict.gates}
        self.assertEqual(by_gate["catalog_crossmatch"], proofs.FAIL)
        self.assertTrue(verdict.is_killed)

    def test_skeptic_controls_run(self):
        self.assertTrue(any("time-shuffle control" in n
                            for n in self.result.skeptic_notes))


class TestLinkingBasics(unittest.TestCase):
    def test_two_detections_cannot_form_a_track(self):
        dets = [
            Detection("a", "epoch_1", 60800.0, 10.0, 20.0, 10, 100),
            Detection("b", "epoch_2", 60801.0, 10.001, 20.0, 10, 100),
        ]
        self.assertEqual(link_tracks(dets), [])

    def test_fit_recovers_rate(self):
        dets = []
        for i, mjd in enumerate([60800.0, 60801.0, 60802.0]):
            dets.append(Detection(f"s{i}", f"epoch_{i + 1}", mjd,
                                  10.0 + 0.001 * i, 20.0, 10, 100))
        fit = fit_motion(dets)
        self.assertAlmostEqual(fit.rate_arcsec_per_day, 3.6 *
                               (1 / 1.0), delta=0.6)


if __name__ == "__main__":
    unittest.main()
