"""Linking, motion-fit, registration, and crew-report evals."""
import math
import unittest

from spherex_hunt import agents
from spherex_hunt.detections import (Detection, angular_sep_arcsec, epoch_offset,
                                     fit_motion, link_tracks, register_epochs)


def det(sid, ep, mjd, ra, dec, snr=25.0):
    return Detection(sid, ep, mjd, ra, dec, snr, 900.0, 0.08, 12.0,
                     0.0, 0.0, "img")


class TestLinking(unittest.TestCase):
    def test_single_epoch_yields_zero_tracks(self):
        dets = [det(f"s{i}", "epoch_1", 60800.0, 10.0 + 0.001 * i, 20.0)
                for i in range(4)]
        self.assertEqual(link_tracks(dets), [])

    def test_tolerance_boundary(self):
        cosd = math.cos(math.radians(20.0))
        step = 5.0 / 3600.0 / cosd  # 5 arcsec/day in RA
        line = [det("a", "epoch_1", 60800.0, 10.0, 20.0),
                det("b", "epoch_2", 60801.0, 10.0 + step, 20.0)]
        on = line + [det("c", "epoch_3", 60802.0, 10.0 + 2 * step, 20.0)]
        off = line + [det("c", "epoch_3", 60802.0, 10.0 + 2 * step,
                          20.0 + 20.0 / 3600.0)]
        self.assertEqual(len(link_tracks(on)), 1)
        self.assertEqual(link_tracks(off), [])

    def test_seed_cap_respected(self):
        pair = [det("a", "epoch_1", 60800.0, 10.0, 20.0),
                det("b", "epoch_2", 60801.0, 10.0 + 2000.0 / 3600.0, 20.0)]
        self.assertEqual(link_tracks(pair, min_epochs=2, max_seed_rate=100.0), [])
        self.assertEqual(len(link_tracks(pair, min_epochs=2)), 1)

    def test_fit_exact_recovery(self):
        dets = [det(f"s{i}", f"epoch_{i + 1}", 60800.0 + i, 10.0 + 0.002 * i, 20.0)
                for i in range(3)]
        fit = fit_motion(dets)
        self.assertAlmostEqual(fit.rate_arcsec_per_day,
                               0.002 * 3600.0 * math.cos(math.radians(20.0)),
                               places=6)
        self.assertAlmostEqual(fit.rms_arcsec, 0.0, places=9)

    def test_angular_sep_values(self):
        self.assertAlmostEqual(angular_sep_arcsec(0.0, 0.0, 1.0, 0.0), 3600.0)
        self.assertAlmostEqual(angular_sep_arcsec(10.0, 20.0, 10.0, 20.0), 0.0)
        # near-pole small offset agrees with the great-circle distance
        got = angular_sep_arcsec(0.0, 89.0, 1.0, 89.0)
        want = 3600.0 * math.degrees(2 * math.asin(
            math.cos(math.radians(89.0)) * math.sin(math.radians(0.5))))
        self.assertAlmostEqual(got, want, delta=1.0)

    def test_epoch_offset_no_match_zeros(self):
        ref = [det("r", "epoch_1", 60800.0, 10.0, 20.0, snr=50.0)]
        mov = [det("m", "epoch_2", 60801.0, 15.0, 25.0, snr=50.0)]
        self.assertEqual(epoch_offset(ref, mov), (0.0, 0.0, 0))

    def test_register_single_epoch_untouched(self):
        dets = [det(f"s{i}", "epoch_1", 60800.0, 10.0 + 0.001 * i, 20.0)
                for i in range(3)]
        out, rep = register_epochs(dets)
        self.assertEqual(out, dets)
        self.assertEqual(rep, {})

    def test_register_flag_note(self):
        dets = [det(f"s{i}", f"epoch_{(i % 3) + 1}", 60800.0 + (i % 3),
                    10.0 + 0.001 * (i % 3), 20.0) for i in range(6)]
        patch = agents.Patch("p", 10.0, 20.0, 0.25)
        reg = agents.run_hunt(patch, dets, register=True)
        self.assertTrue(any("registration offsets" in n for n in reg.skeptic_notes))
        plain = agents.run_hunt(patch, dets)
        self.assertFalse(any("registration offsets" in n for n in plain.skeptic_notes))

    def test_skeptic_deterministic(self):
        dets = [det(f"s{i}", f"epoch_{(i % 4) + 1}", 60800.0 + (i % 4),
                    10.0 + 0.001 * (i % 4), 20.0 + 0.0005 * (i % 4))
                for i in range(8)]
        patch = agents.Patch("p", 10.0, 20.0, 0.25)
        a = agents.run_hunt(patch, dets).skeptic_notes
        b = agents.run_hunt(patch, dets).skeptic_notes
        self.assertEqual(a, b)
        self.assertTrue(any("time-shuffle control" in n for n in a))
        self.assertTrue(any("position-scramble control" in n for n in a))

    def test_reporter_and_case_file(self):
        dets = [det("m", f"epoch_{i + 1}", 60800.0 + i, 10.0 + 0.001 * i, 20.0,
                    snr=50.0) for i in range(3)]
        patch = agents.Patch("p", 10.0, 20.0, 0.25)
        res = agents.run_hunt(patch, dets, catalog=[(100.0, -20.0, "far")])
        text = agents.reporter(res)
        self.assertIn("Tracks linked:", text)
        self.assertIn("Case files for human review:", text)
        md = agents.case_file_markdown(res, 0)
        self.assertIn("# Case file:", md)
        self.assertIn("## Proof gates", md)
        self.assertIn("## Detections (provenance)", md)


if __name__ == "__main__":
    unittest.main()
