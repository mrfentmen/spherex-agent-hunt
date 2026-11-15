import unittest

from spherex_hunt.sia import ImageRecord, detector_key, pick_overlapping_frames


class TestOverlap(unittest.TestCase):
    def test_picks_same_sky_spaced_epochs(self):
        def rec(oid, ra, t):
            return ImageRecord(oid, "spherex_qr3", ra, 54.35, t, t + 0.001,
                               f"https://example.invalid/{oid}")

        recs = [
            rec("a", 210.80, 61177.0),
            rec("a2", 211.60, 61177.0),
            rec("b", 210.801, 61178.2),
            rec("b2", 210.800, 61178.3),
            rec("c", 210.799, 61180.0),
        ]
        picks = pick_overlapping_frames(recs, min_gap_days=1.0, max_sep_deg=0.5, limit=3)
        self.assertEqual([r.obs_id for r in picks], ["a", "b", "c"])

    def test_picks_single_detector(self):
        def rec(oid, ra, t, det):
            return ImageRecord(oid, "spherex_qr3", ra, 54.35, t, t + 0.001,
                               f"https://irsa.ipac.caltech.edu/ibe/data/x/level2_{oid}_{det}_spx.fits")

        recs = [
            rec("a", 210.80, 61177.0, "1D2"),
            rec("b", 210.801, 61178.2, "1D5"),
            rec("c", 210.799, 61180.0, "1D2"),
        ]
        picks = pick_overlapping_frames(recs, min_gap_days=1.0, max_sep_deg=0.5,
                                        limit=3, same_detector=True)
        self.assertEqual([r.obs_id for r in picks], ["a", "c"])
        self.assertEqual(detector_key(picks[0].access_url), "1D2")


if __name__ == "__main__":
    unittest.main()
