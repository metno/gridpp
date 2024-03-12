from __future__ import print_function
import unittest
import gridpp
import numpy as np
import os


class Test(unittest.TestCase):
    def test_invalid_input(self):
        """Check that dimension missmatch results in error"""
        with self.assertRaises(Exception) as e:
            gridpp.qnh([101325], [0, 20])

    def test_invalid_values(self):
        self.assertTrue(np.isnan(gridpp.qnh([-1], [0])))
        self.assertTrue(np.isnan(gridpp.qnh([101325], [np.nan])))
        self.assertTrue(np.isnan(gridpp.qnh([np.nan], [0])))

    def test_1(self):
        p = [101325, 90000, 90000, 110000]
        alt = [0, 1000, 0, -1000]
        expected = [101325, 101463.21875, 90000, 97752.90742927508]
        for i in range(len(p)):
            self.assertAlmostEqual(gridpp.qnh(p[i], alt[i]), expected[i], 1)
        np.testing.assert_almost_equal(gridpp.qnh(p, alt), expected, 1)

    def test_no_pressure(self):
        for altitude in [-1000, 0, 1000]:
            self.assertEqual(gridpp.qnh([0], [altitude]), [0])
            self.assertEqual(gridpp.qnh(0, altitude), 0)

    def test_empty(self):
        np.testing.assert_almost_equal(gridpp.qnh([],[]), [])


if __name__ == '__main__':
    unittest.main()
