import unittest

import numpy as np

from spherex_hunt.extract import find_peaks, local_snr, read_image_and_variance


class TestExtract(unittest.TestCase):
    def test_find_peaks_two_blobs(self):
        img = np.zeros((64, 64), dtype=np.float32)
        img[10, 10] = 100.0
        img[40, 30] = 90.0
        peaks = find_peaks(img, nsigma=5.0)
        self.assertEqual(len(peaks), 2)

    def test_plateau_yields_one_peak(self):
        img = np.zeros((32, 32), dtype=np.float32)
        img[10:13, 10:13] = 50.0
        peaks = find_peaks(img, nsigma=1.0)
        self.assertEqual(len(peaks), 1)

    def test_edge_peaks_rejected(self):
        img = np.zeros((32, 32), dtype=np.float32)
        img[1, 1] = 100.0
        img[30, 30] = 90.0
        peaks = find_peaks(img, nsigma=1.0, edge=5)
        self.assertEqual(len(peaks), 0)

    def test_local_snr_prefers_variance(self):
        self.assertAlmostEqual(local_snr(10.0, 2.0, 4.0, 99.0), 4.0)
        self.assertAlmostEqual(local_snr(10.0, 2.0, 0.0, 4.0), 2.0)
        self.assertAlmostEqual(local_snr(10.0, 2.0, float("nan"), 4.0), 2.0)

    def test_band_from_path(self):
        from spherex_hunt.extract import band_from_path
        self.assertEqual(band_from_path("/x/level2_2026W20_2A_0464_1D2_spx.fits"), "1D2")
        self.assertEqual(band_from_path("nodetector.fits"), "")

    def test_aperture_flux_subtracts_background(self):
        from spherex_hunt.extract import aperture_flux
        img = np.full((32, 32), 10.0, dtype=np.float32)
        img[16, 16] = 110.0
        self.assertAlmostEqual(aperture_flux(img, 16.0, 16.0), 100.0, delta=15.0)

    def test_nan_region_peak_rejected(self):
        img = np.zeros((32, 32), dtype=np.float32)
        img[5:25, 5:25] = np.nan
        img[15, 15] = 50.0
        img[28, 28] = 50.0
        peaks = find_peaks(img, nsigma=1.0, edge=2)
        self.assertEqual([(x, y) for x, y, _ in peaks], [(28, 28)])

    def test_read_image_and_variance_roundtrip(self):
        import os
        import tempfile

        def card(k, v):
            return (k.ljust(8) + "= " + v.rjust(20)).ljust(80)[:80]

        def hdu(pairs, data=b""):
            blk = "".join([card(k, v) for k, v in pairs] + ["END".ljust(80)])
            while len(blk) % 2880:
                blk += " " * 80
            raw = blk.encode("ascii") + data
            return raw + b"\0" * ((-len(raw)) % 2880)

        img = np.arange(16, dtype=np.float32).reshape(4, 4)
        var = np.full((4, 4), 2.0, dtype=np.float32)
        primary = hdu([("SIMPLE", "T"), ("BITPIX", "8"),
                       ("NAXIS", "0"), ("EXTEND", "T")])
        image = hdu([("XTENSION", "'IMAGE   '"), ("BITPIX", "-32"),
                     ("NAXIS", "2"), ("NAXIS1", "4"),
                     ("NAXIS2", "4"), ("EXTNAME", "'IMAGE   '")],
                    img.astype(">f4").tobytes())
        variance = hdu([("XTENSION", "'IMAGE   '"), ("BITPIX", "-32"),
                        ("NAXIS", "2"), ("NAXIS1", "4"),
                        ("NAXIS2", "4"), ("EXTNAME", "'VARIANCE'")],
                       var.astype(">f4").tobytes())
        with tempfile.TemporaryDirectory() as tmp:
            path = os.path.join(tmp, "t.fits")
            with open(path, "wb") as fh:
                fh.write(primary + image + variance)
            hdr, got, got_var = read_image_and_variance(path)
        self.assertEqual(got.shape, (4, 4))
        self.assertTrue(np.allclose(got, img))
        self.assertTrue(np.allclose(got_var, var))


if __name__ == '__main__':
    unittest.main()
