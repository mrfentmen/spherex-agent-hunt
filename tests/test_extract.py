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


if __name__ == '__main__':
    unittest.main()
