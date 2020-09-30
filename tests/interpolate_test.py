from __future__ import print_function
import unittest
import gridpp
import numpy as np


class Test(unittest.TestCase):
    def test_basic(self):
        x = [0, 1, 2]
        y = [0, 2, 1]

        # Test numbers in interpolator points
        self.assertEqual(gridpp.interpolate(0, x, y), 0)
        self.assertEqual(gridpp.interpolate(2, x, y), 1)
        self.assertEqual(gridpp.interpolate(1, x, y), 2)

        # Test numbers in between interpolator points
        self.assertAlmostEqual(gridpp.interpolate(0.5, x, y), 1)
        self.assertAlmostEqual(gridpp.interpolate(0.9, x, y), 1.8)
        self.assertAlmostEqual(gridpp.interpolate(1.5, x, y), 1.5)

        # Test outside range
        self.assertEqual(gridpp.interpolate(-1, x, y), 0)
        self.assertEqual(gridpp.interpolate(3, x, y), 1)

    def test_single(self):
        x = [0]
        y = [0]

        self.assertEqual(gridpp.interpolate(0, x, y), 0)
        self.assertEqual(gridpp.interpolate(-1, x, y), 0)
        self.assertEqual(gridpp.interpolate(1, x, y), 0)

    def test_empty(self):
        x = []
        y = []

        self.assertTrue(np.isnan(gridpp.interpolate(0, x, y)))
        self.assertTrue(np.isnan(gridpp.interpolate(-1, x, y)))
        self.assertTrue(np.isnan(gridpp.interpolate(1, x, y)))

if __name__ == '__main__':
    unittest.main()
