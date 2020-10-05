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
        quantiles = [-0.1, 1.1, np.nan, np.inf]

        for quantile in quantiles:
            with self.assertRaises(ValueError) as e:
                gridpp.neighbourhood_quantile(field, quantile, radius)

    def test_empty(self):
        for quantile in np.arange(0.1,0.9,0.1):
            for num_thresholds in [1, 2]:
                thresholds = gridpp.get_neighbourhood_thresholds(values, num_thresholds)
                output = gridpp.neighbourhood_quantile_fast([[]], 0.9, 1, thresholds)
                self.assertEqual(len(output.shape), 2)
                self.assertEqual(output.shape[0], 0)
                self.assertEqual(output.shape[1], 0)

    def test_missing(self):
        empty = np.zeros([5, 5])
        empty[0:3, 0:3] = np.nan
        output = gridpp.neighbourhood_quantile(empty, 0.5, 1)
        self.assertTrue(np.isnan(np.array(output)[0:2,0:2]).all())

    def test_quantile(self):
        thresholds = gridpp.get_neighbourhood_thresholds(values, 100)
        output = np.array(gridpp.neighbourhood_quantile_fast(values, 0.5, 1, thresholds))
        self.assertEqual(output[2][2], 12)   # Should be 12.5
        self.assertEqual(output[2][3], 12.5) # Should be 13

        output = np.array(gridpp.neighbourhood_quantile_fast(np.full([100,100], np.nan), 0.5, 1, thresholds))
        self.assertTrue(np.isnan(np.array(output)).all())

        output = np.array(gridpp.neighbourhood_quantile_fast(np.zeros([100,100]), 0.5, 1, thresholds))
        self.assertTrue((np.array(output) == 0).all())

        output = np.array(gridpp.neighbourhood_quantile(values, 0.5, 1))
        self.assertEqual(output[2][2], 12.5)
        self.assertEqual(output[2][3], 13)
        self.assertEqual(output[0][4], 4)

    def test_quantile_exceed(self):
        input = np.reshape(range(25), [5, 5]).astype(float)
        thresholds = [1, 5, 10]
        output = gridpp.neighbourhood_quantile_fast(values, 0.9, 1, thresholds)

    def test_mean(self):
        input = (np.random.rand(1000, 1000)> 0.5).astype(float)
        thresholds = [1, 5, 10]
        output = gridpp.neighbourhood(values, 7, gridpp.Mean)

if __name__ == '__main__':
    unittest.main()
