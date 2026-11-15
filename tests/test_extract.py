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

    def test_sip_linear_term_matches_cd(self):
        from spherex_hunt.extract import sip_to_radec
        hdr = {"CRVAL1": "10.0", "CRVAL2": "20.0", "CRPIX1": "0.0", "CRPIX2": "0.0",
               "CD1_1": "0.001", "CD1_2": "0.0", "CD2_1": "0.0", "CD2_2": "0.001",
               "A_ORDER": "0", "B_ORDER": "0"}
        ra, dec = sip_to_radec(hdr, 100.0, 0.0)
        import math
        self.assertAlmostEqual(ra, 10.0 + 0.1 / math.cos(math.radians(20.0)), places=5)
        self.assertAlmostEqual(dec, 20.0, places=4)  # real TAN, not flat-sky

    def test_sip_quadratic_term_applies(self):
        from spherex_hunt.extract import sip_to_radec
        hdr = {"CRVAL1": "0.0", "CRVAL2": "0.0", "CRPIX1": "0.0", "CRPIX2": "0.0",
               "CD1_1": "1.0", "CD1_2": "0.0", "CD2_1": "0.0", "CD2_2": "1.0",
               "A_ORDER": "2", "B_ORDER": "0", "A_2_0": "0.001"}
        ra, dec = sip_to_radec(hdr, 0.1, 0.0)
        self.assertAlmostEqual(ra, 0.1 + 0.001 * 0.01, places=6)
        self.assertAlmostEqual(dec, 0.0, places=6)

    def test_sip_constant_term_applies(self):
        from spherex_hunt.extract import sip_to_radec
        hdr = {"CRVAL1": "0.0", "CRVAL2": "0.0", "CRPIX1": "0.0", "CRPIX2": "0.0",
               "CD1_1": "1.0", "CD1_2": "0.0", "CD2_1": "0.0", "CD2_2": "1.0",
               "A_ORDER": "1", "B_ORDER": "0", "A_0_0": "0.5"}
        ra, dec = sip_to_radec(hdr, 0.0, 0.0)
        self.assertAlmostEqual(ra, 0.5, places=4)  # exact TAN, not flat-sky
        self.assertAlmostEqual(dec, 0.0, places=6)

    def test_sip_reverse_roundtrip(self):
        from spherex_hunt.extract import sip_reverse, sip_to_radec
        hdr = {"CRVAL1": "10.0", "CRVAL2": "20.0", "CRPIX1": "100.0", "CRPIX2": "100.0",
               "CD1_1": "0.001", "CD1_2": "0.0002", "CD2_1": "-0.0001", "CD2_2": "0.0011",
               "A_ORDER": "2", "B_ORDER": "2", "A_2_0": "0.001", "B_0_2": "-0.002",
               "AP_ORDER": "2", "BP_ORDER": "2"}
        # reverse of the forward map needs the AP/BP pair; here only the
        # linear part round-trips exactly, distortion to first order
        ra, dec = sip_to_radec(hdr, 105.0, 95.0)
        x, y = sip_reverse(hdr, ra, dec)
        self.assertAlmostEqual(x, 105.0, delta=2.0)
        self.assertAlmostEqual(y, 95.0, delta=2.0)

    def test_all_nan_yields_no_peaks(self):
        img = np.full((32, 32), np.nan, dtype=np.float32)
        self.assertEqual(find_peaks(img), [])

    def test_flat_image_yields_no_peaks(self):
        img = np.full((32, 32), 7.0, dtype=np.float32)
        self.assertEqual(find_peaks(img), [])

    def test_aperture_nan_and_border(self):
        from spherex_hunt.extract import aperture_flux
        img = np.full((32, 32), np.nan, dtype=np.float32)
        self.assertEqual(aperture_flux(img, 16.0, 16.0), 0.0)
        img2 = np.full((32, 32), 10.0, dtype=np.float32)
        img2[0, 0] = 110.0
        self.assertGreater(aperture_flux(img2, 0.0, 0.0), 0.0)

    def test_moments_subpixel_centroid(self):
        from spherex_hunt.extract import moments
        yy, xx = np.mgrid[0:32, 0:32]
        img = (1000.0 * np.exp(-((xx - 10.3) ** 2 + (yy - 20.7) ** 2) / 2.0)).astype(np.float32)
        cx, cy, _, _, _ = moments(img, 10, 21)
        self.assertAlmostEqual(cx, 10.3, delta=0.1)
        self.assertAlmostEqual(cy, 20.7, delta=0.1)

    def test_variance_absent_returns_none(self):
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

        img = np.ones((4, 4), dtype=np.float32)
        primary = hdu([("SIMPLE", "T"), ("BITPIX", "8"),
                       ("NAXIS", "0"), ("EXTEND", "T")])
        image = hdu([("XTENSION", "'IMAGE   '"), ("BITPIX", "-32"),
                     ("NAXIS", "2"), ("NAXIS1", "4"),
                     ("NAXIS2", "4"), ("EXTNAME", "'IMAGE   '")],
                    img.astype(">f4").tobytes())
        with tempfile.TemporaryDirectory() as tmp:
            path = os.path.join(tmp, "t.fits")
            with open(path, "wb") as fh:
                fh.write(primary + image)
            _, got, got_var = read_image_and_variance(path)
        self.assertEqual(got.shape, (4, 4))
        self.assertIsNone(got_var)


if __name__ == '__main__':
    unittest.main()
