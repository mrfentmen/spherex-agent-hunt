"""Locked ground-truth receipts: calibration numbers must never drift."""
import os
import sys
import tempfile
import unittest

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "src"))

from spherex_hunt import agents, calibration, proofs  # noqa: E402


class TestReceipts(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.tmp = tempfile.mkdtemp()
        cls.csv_path = os.path.join(cls.tmp, "calibration.csv")
        calibration.generate(cls.csv_path)
        cls.detections = calibration.load_calibration(cls.csv_path)
        patch = agents.Patch("receipt-patch", 210.80225, 54.34894, 0.25)
        cls.result = agents.run_hunt(patch, cls.detections)

    def mover_track(self):
        for track, verdict, _ in self.result.tracks:
            if any(d.source_id == "PLANTED-MOVER" for d in track.detections):
                return track, verdict
        self.fail("planted mover track missing")

    def test_receipt_numbers(self):
        track, _ = self.mover_track()
        self.assertAlmostEqual(track.fit.rate_arcsec_per_day, 32.33, delta=0.5)
        self.assertEqual(track.epochs, 4)
        self.assertLess(track.fit.rms_arcsec, 0.2)

    def test_receipt_holdout(self):
        _, verdict = self.mover_track()
        holdout = {g.gate: g for g in verdict.gates}["holdout_prediction"]
        self.assertEqual(holdout.status, proofs.PASS)
        self.assertIn("0.26", holdout.detail)

    def test_receipt_zero_case_files(self):
        self.assertEqual(sum(1 for _, v, _ in self.result.tracks if v.is_case_file), 0)

    def test_stars_killed_as_stationary(self):
        star_tracks = [t for t, _, _ in self.result.tracks
                       if any(d.source_id.startswith("STAR-") for d in t.detections)]
        self.assertGreaterEqual(len(star_tracks), 6)
        for track in star_tracks:
            self.assertLess(track.fit.rate_arcsec_per_day, 1.0)

    def test_glitches_killed_on_artifact(self):
        for track, verdict, _ in self.result.tracks:
            ids = {d.source_id for d in track.detections}
            if ids & {"PLANTED-MOVER"} or any(i.startswith("STAR-") for i in ids):
                continue
            if any(i.startswith("GLITCH-") for i in ids):
                failed = [g.gate for g in verdict.gates if g.status == proofs.FAIL]
                self.assertIn("artifact_exclusion", failed)


if __name__ == "__main__":
    unittest.main()
