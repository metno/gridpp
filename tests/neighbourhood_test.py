from __future__ import print_function
import unittest
import gridpp
import numpy as np


lats = [60, 60, 60, 60, 60, 70]
lons = [10,10.1,10.2,10.3,10.4, 10]

"""Simple check
20 21 22 23 24
15 16 17 18 19
10 11 12 13 nan
5  6  7  nan  9
0  1  2  3  4
"""
values = np.reshape(range(25), [5, 5]).astype(float)
values[1, 3] = np.nan
values[2, 4] = np.nan
values = np.array(values)

class Test(unittest.TestCase):
    def test_invalid_arguments(self):
        """Check that exception is thrown for invalid arguments"""
        field = np.ones([5, 5])
        radius = -1
        stats = [gridpp.Mean, gridpp.Min, gridpp.Max, gridpp.Median]

        for stat in stats:
            with self.assertRaises(ValueError) as e:
                gridpp.neighbourhood(field, radius, stat)

        # User should use the _quantile function
        with self.assertRaises(Exception) as e:
            gridpp.neighbourhood(field, 1, gridpp.Quantile)

    def test_empty(self):
        """Empty input array"""
        for statistic in [gridpp.Mean, gridpp.Min, gridpp.Max, gridpp.Median, gridpp.Std, gridpp.Variance]:
            output = gridpp.neighbourhood([[]], 1, statistic)
            self.assertEqual(len(output.shape), 2)
            self.assertEqual(output.shape[0], 0)
            self.assertEqual(output.shape[1], 0)

    def test_missing(self):
        """Missing values in input array"""
        empty = np.zeros([5, 5])
        empty[0:3, 0:3] = np.nan
        for statistic in [gridpp.Mean, gridpp.Min, gridpp.Max, gridpp.Median, gridpp.Std, gridpp.Variance]:
            output = gridpp.neighbourhood(empty, 1, statistic)
            self.assertTrue(np.isnan(np.array(output)[0:2,0:2]).all())

    def test_mean(self):
        output = gridpp.neighbourhood(values, 1, gridpp.Mean)
        self.assertEqual(output[2][2], 12.5)
        self.assertAlmostEqual(output[0][4], 5.3333, 4)

        output = gridpp.neighbourhood(values, 100, gridpp.Mean)
        self.assertTrue((np.abs(np.array(output) - 12.086956)<0.0001).all())

        output = np.array(gridpp.neighbourhood(values, 0, gridpp.Mean)).flatten()
        I = np.where(np.isnan(output) == 0)[0]
        self.assertTrue((np.isnan(output) == np.isnan(values.flatten())).all())
        self.assertTrue((output[I] == values.flatten()[I]).all())

    def test_min(self):
        output = gridpp.neighbourhood(values, 1, gridpp.Min)
        self.assertEqual(output[2][2], 6)
        output = gridpp.neighbourhood(values, 1, gridpp.Min)
        self.assertEqual(output[0][4], 3)

        output = gridpp.neighbourhood(values, 100, gridpp.Min)
        self.assertTrue((np.array(output) == 0).all())

    def test_max(self):
        output = gridpp.neighbourhood(values, 1, gridpp.Max)
        self.assertEqual(output[2][2], 18)
        output = gridpp.neighbourhood(values, 1, gridpp.Max)
        self.assertEqual(output[0][4], 9)

        output = gridpp.neighbourhood(values, 100, gridpp.Max)
        self.assertTrue((np.array(output) == 24).all())

    def test_mean(self):
        input = (np.random.rand(1000, 1000)> 0.5).astype(float)
        thresholds = [1, 5, 10]
        output = gridpp.neighbourhood(values, 7, gridpp.Mean)

if __name__ == '__main__':
    unittest.main()
