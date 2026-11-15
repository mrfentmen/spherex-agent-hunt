import unittest

from spherex_hunt.catalogs import parse_tap_csv

SAMPLE = """ra,dec,phot_g_mean_mag,source_id
210.801,54.351,15.2,123456789
210.802,54.352,18.7,987654321
"""


class TestCatalogs(unittest.TestCase):
    def test_parse_tap_csv(self):
        rows = parse_tap_csv(SAMPLE)
        self.assertEqual(len(rows), 2)
        self.assertAlmostEqual(rows[0][0], 210.801)
        self.assertAlmostEqual(rows[0][1], 54.351)
        self.assertIn("123456789", rows[0][2])

    def test_parse_empty(self):
        self.assertEqual(parse_tap_csv("ra,dec\n"), [])


if __name__ == "__main__":
    unittest.main()
