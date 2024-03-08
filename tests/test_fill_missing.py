from __future__ import print_function
import unittest
import gridpp
import numpy as np


class Test(unittest.TestCase):
    def test_linear(self):
        """Check that we are able to recover the missing values"""
        values0 = np.reshape(np.arange(25), [5, 5]).astype(float)

        # Add some missing values
        values = np.copy(values0)
        values[2, 1:4] = np.nan
        values[1, 1] = np.nan
        output = gridpp.fill_missing(values)
        np.testing.assert_array_equal(output, values0)

    def test_missing_on_edge(self):
        """Check that we can recover when one dimension has missing all the way to the edge"""
        values0 = np.reshape(np.arange(25), [5, 5]).astype(float)

        # Add some missing values
        values = np.copy(values0)
        values[1, 1] = np.nan
        values[1, 3:5] = np.nan
        values[1, 4] = np.nan
        values[1, 0:2] = np.nan
        output = gridpp.fill_missing(values)
        np.testing.assert_array_equal(output, values0)

    def test_missing_on_y_edge(self):
        """Regression test for bug when X is wider than Y, and a y-slice has missing all the way to the upper edge"""
        values0 = np.reshape(np.arange(24), [3, 8]).astype(float)
        values = np.copy(values0)
        values[1:, 1] = np.nan
        output = gridpp.fill_missing(values)
        np.testing.assert_array_equal(output, values0)

    def test_missing_on_both_edges(self):
        """Check that we are able to recover the missing values"""
        values0 = np.reshape(np.arange(25), [5, 5]).astype(float)

        # Add some missing values
        values = np.copy(values0)
        values[3:5, 3:5] = np.nan
        output = gridpp.fill_missing(values)
        np.testing.assert_array_equal(output[0:3, :], values0[0:3, :])
        np.testing.assert_array_equal(output[:, 0:3], values0[:, 0:3])
        np.testing.assert_array_equal(output[3:5, 3:5], np.nan * np.zeros([2, 2]))


if __name__ == '__main__':
    unittest.main()
