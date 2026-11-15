"""Gate-by-gate proof evals: every PASS/FAIL/UNVERIFIED edge."""
import unittest

from spherex_hunt import proofs
from spherex_hunt.detections import Detection, Track, fit_motion


def det(sid, ep, mjd, ra, dec, snr=25.0, flux=900.0, rnd=0.08, area=12.0,
        band="", fb=None):
    return Detection(sid, ep, mjd, ra, dec, snr, flux, rnd, area, 0.0, 0.0,
                     "img", band, fb)


def mover_dets(n=3, dra=0.001, dec0=20.0, **kw):
    return [det("m", f"epoch_{i + 1}", 60800.0 + i, 10.0 + dra * i, dec0, **kw)
            for i in range(n)]


def track_of(dets):
    return Track(detections=dets, fit=fit_motion(dets))


def gate(track, name, **kw):
    v = proofs.evaluate(track, **kw)
    return {g.gate: g.status for g in v.gates}[name]


class TestGates(unittest.TestCase):
    def test_epochs_two_fail_three_pass(self):
        self.assertEqual(gate(track_of(mover_dets(2)), "repeat_epochs"), proofs.FAIL)
        self.assertEqual(gate(track_of(mover_dets(3)), "repeat_epochs"), proofs.PASS)

    def test_linear_kills_stationary_and_absurd(self):
        still = [det("s", f"epoch_{i + 1}", 60800.0 + i, 10.0, 20.0) for i in range(3)]
        self.assertEqual(gate(track_of(still), "linear_motion"), proofs.FAIL)
        fast = [det("f", f"epoch_{i + 1}", 60800.0 + i, 10.0 + 0.6 * i, 20.0)
                for i in range(3)]
        self.assertEqual(gate(track_of(fast), "linear_motion"), proofs.FAIL)
        self.assertEqual(gate(track_of(mover_dets(3)), "linear_motion"), proofs.PASS)

    def test_linear_rms_boundary(self):
        ok = [det("m", f"epoch_{i + 1}", 60800.0 + i, 10.0 + 0.001 * i,
                  20.0 + (0.0007 if i % 2 else -0.0007)) for i in range(4)]
        self.assertEqual(gate(track_of(ok), "linear_motion"), proofs.PASS)
        bad = [det("m", f"epoch_{i + 1}", 60800.0 + i, 10.0 + 0.001 * i,
                   20.0 + (0.001 if i % 2 else -0.001)) for i in range(4)]
        self.assertEqual(gate(track_of(bad), "linear_motion"), proofs.FAIL)

    def test_holdout_paths(self):
        self.assertEqual(gate(track_of(mover_dets(3)), "holdout_prediction"), proofs.PASS)
        off = mover_dets(3)
        off[-1] = det("m", "epoch_3", 60802.0, 10.0 + 0.003 + 0.01, 20.0)
        self.assertEqual(gate(track_of(off), "holdout_prediction"), proofs.FAIL)
        self.assertEqual(gate(track_of(mover_dets(2)), "holdout_prediction"),
                         proofs.UNVERIFIED)

    def test_artifact_each_kill(self):
        self.assertEqual(gate(track_of(mover_dets(3)), "artifact_exclusion"), proofs.PASS)
        self.assertEqual(gate(track_of(mover_dets(3, rnd=0.9)), "artifact_exclusion"),
                         proofs.FAIL)
        self.assertEqual(gate(track_of(mover_dets(3, area=1.0)), "artifact_exclusion"),
                         proofs.FAIL)
        self.assertEqual(gate(track_of(mover_dets(3, snr=1.0)), "artifact_exclusion"),
                         proofs.FAIL)

    def test_catalog_paths(self):
        tr = track_of(mover_dets(3))
        self.assertEqual(gate(tr, "catalog_crossmatch"), proofs.UNVERIFIED)
        self.assertEqual(gate(tr, "catalog_crossmatch",
                              catalog=[(100.0, -20.0, "far")]), proofs.PASS)
        fit = tr.fit
        self.assertEqual(gate(tr, "catalog_crossmatch",
                              catalog=[(fit.ra0, fit.dec0, "known")]), proofs.FAIL)

    def test_spectrum_paths(self):
        tr = track_of(mover_dets(3))
        self.assertEqual(gate(tr, "spectrum_consistency"), proofs.UNVERIFIED)
        steady = [det("m", f"epoch_{i + 1}", 60800.0 + i, 10.0 + 0.001 * i, 20.0,
                      band="D2", fb=900.0 + i) for i in range(3)]
        self.assertEqual(gate(track_of(steady), "spectrum_consistency"), proofs.PASS)
        wild = [det("m", f"epoch_{i + 1}", 60800.0 + i, 10.0 + 0.001 * i, 20.0,
                    band="D2", fb=v) for i, v in enumerate([100.0, 900.0, 100.0])]
        self.assertEqual(gate(track_of(wild), "spectrum_consistency"), proofs.FAIL)

    def test_verdict_labels(self):
        steady = [det("m", f"epoch_{i + 1}", 60800.0 + i, 10.0 + 0.001 * i, 20.0,
                      band="D2", fb=900.0) for i in range(3)]
        v = proofs.evaluate(track_of(steady), catalog=[(100.0, -20.0, "far")])
        self.assertTrue(v.is_case_file)
        self.assertTrue(v.label.startswith("CASE FILE"))
        still = [det("s", f"epoch_{i + 1}", 60800.0 + i, 10.0, 20.0) for i in range(3)]
        vk = proofs.evaluate(track_of(still))
        self.assertTrue(vk.is_killed)
        self.assertTrue(vk.label.startswith("KILLED"))
        vc = proofs.evaluate(track_of(mover_dets(3)))
        self.assertFalse(vc.is_case_file)
        self.assertFalse(vc.is_killed)
        self.assertIn("CANDIDATE", vc.label)


if __name__ == "__main__":
    unittest.main()
