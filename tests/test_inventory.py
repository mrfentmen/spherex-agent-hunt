"""Inventory, picker, CSV, and empty-hunt evals."""
import csv
import os
import tempfile
import unittest

from spherex_hunt import agents, sia
from spherex_hunt.detections import Detection, load_detections


def rec(oid, ra, t, url="https://example.invalid/x_1D2_spx.fits"):
    return sia.ImageRecord(oid, "spherex_qr3", ra, 54.35, t, t + 0.001, url)


class TestInventory(unittest.TestCase):
    def test_no_table_raises(self):
        bad = b"""<?xml version="1.0"?><VOTABLE xmlns="http://www.ivoa.net/xml/VOTable/v1.3"><RESOURCE></RESOURCE></VOTABLE>"""
        with self.assertRaises(ValueError):
            sia.parse_votable(bad)

    def test_bad_floats_become_none(self):
        doc = b"""<?xml version="1.0"?>
<VOTABLE xmlns="http://www.ivoa.net/xml/VOTable/v1.3"><RESOURCE><TABLE>
<FIELD name="obs_id"/><FIELD name="obs_collection"/><FIELD name="s_ra"/>
<FIELD name="s_dec"/><FIELD name="t_min"/><FIELD name="t_max"/><FIELD name="access_url"/>
<DATA><TABLEDATA>
<TR><TD>img-1</TD><TD>spherex_qr3</TD><TD>not-a-number</TD><TD></TD>
<TD>60800.0</TD><TD>60800.01</TD><TD></TD></TR>
</TABLEDATA></DATA></TABLE></RESOURCE></VOTABLE>"""
        (r,) = sia.parse_votable(doc)
        self.assertIsNone(r.ra_deg)
        self.assertIsNone(r.dec_deg)
        self.assertEqual(r.access_url, "")

    def test_detector_variants(self):
        self.assertEqual(sia.detector_key("a_3D6_spx.fits"), "3D6")
        self.assertEqual(sia.detector_key("a_2D5_spx.fits"), "2D5")
        self.assertEqual(sia.detector_key("plain.fits"), "")
        self.assertEqual(sia.detector_key(""), "")
        # filenames are uppercase tags; lowercase does not match
        self.assertEqual(sia.detector_key("a_1d2_spx.fits"), "")

    def test_picker_empty_and_limit(self):
        self.assertEqual(sia.pick_overlapping_frames([]), [])
        recs = [rec(f"i{i}", 210.80 + 0.001 * i, 61177.0 + 2.0 * i) for i in range(5)]
        picks = sia.pick_overlapping_frames(recs, min_gap_days=1.0, max_sep_deg=0.5, limit=2)
        self.assertEqual(len(picks), 2)

    def test_picker_order_and_gaps(self):
        recs = [rec("late", 210.801, 61180.0),
                rec("mid", 210.800, 61178.0),
                rec("tooclose", 210.799, 61178.1),
                rec("early", 210.798, 61177.0)]
        picks = sia.pick_overlapping_frames(recs, min_gap_days=1.0, max_sep_deg=0.5, limit=3)
        self.assertEqual([r.obs_id for r in picks], ["early", "mid", "late"])

    def test_hunt_empty(self):
        patch = agents.Patch("empty", 10.0, 20.0, 0.25)
        res = agents.run_hunt(patch, [])
        self.assertEqual(res.tracks, [])
        self.assertIn("Tracks linked: 0", agents.reporter(res))

    def test_csv_round_trip(self):
        dets = [Detection("s1", "epoch_1", 60800.0, 10.5, 20.5, 30.0, 500.0,
                          0.1, 9.0, 100.0, 200.0, "obs-1", "1D2", 480.0)]
        cols = ["source_id", "epoch_id", "mjd", "ra_deg", "dec_deg", "snr",
                "flux", "roundness", "pixel_area", "x_px", "y_px", "image_id",
                "band", "flux_band"]
        with tempfile.TemporaryDirectory() as tmp:
            path = os.path.join(tmp, "d.csv")
            with open(path, "w", newline="") as fh:
                w = csv.DictWriter(fh, fieldnames=cols)
                w.writeheader()
                for d in dets:
                    w.writerow({c: getattr(d, c) for c in cols})
            back = load_detections(path)
        self.assertEqual(len(back), 1)
        for c in cols:
            self.assertEqual(getattr(back[0], c), getattr(dets[0], c))

    def test_load_tolerates_missing_optionals(self):
        with tempfile.TemporaryDirectory() as tmp:
            path = os.path.join(tmp, "d.csv")
            with open(path, "w", newline="") as fh:
                fh.write("source_id,epoch_id,mjd,ra_deg,dec_deg,snr,flux,roundness,pixel_area,x_px,y_px,image_id,band,flux_band\n")
                fh.write("s1,epoch_1,60800.0,10.0,20.0,10,100,0.1,9,0,0,img,,\n")
            (d,) = load_detections(path)
        self.assertEqual(d.band, "")
        self.assertIsNone(d.flux_band)

    def test_inventory_json_round_trip(self):
        recs = [rec("a", 210.8, 61177.0)]
        with tempfile.TemporaryDirectory() as tmp:
            path = os.path.join(tmp, "inv.json")
            sia.write_inventory(recs, path)
            import json
            (row,) = json.load(open(path))
        self.assertEqual(row["obs_id"], "a")
        self.assertAlmostEqual(row["ra_deg"], 210.8)


if __name__ == "__main__":
    unittest.main()
