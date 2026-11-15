import unittest
import numpy as np
from spherex_hunt.extract import find_peaks


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


if __name__ == '__main__':
    unittest.main()
