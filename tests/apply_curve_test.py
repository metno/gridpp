from __future__ import print_function
import unittest
import gridpp
import time
import numpy as np
import sys

class Test(unittest.TestCase):
    def test_empty_curve(self):
        """Check for exception on empty curve"""
        with self.assertRaises(Exception) as e:
            gridpp.apply_curve([0, 1], [[], []], gridpp.OneToOne, gridpp.OneToOne)
        with self.assertRaises(Exception) as e:
            gridpp.apply_curve([0, 1], [[1, 2], []], gridpp.OneToOne, gridpp.OneToOne)
        with self.assertRaises(Exception) as e:
            gridpp.apply_curve([0, 1], [[], [1, 2]], gridpp.OneToOne, gridpp.OneToOne)

    def test_invalid_curve(self):
        """Check for exception on invalid curve"""
        with self.assertRaises(Exception) as e:
            gridpp.apply_curve([0, 1], [[1, 2, 3], [1, 2]], gridpp.OneToOne, gridpp.OneToOne)
        with self.assertRaises(Exception) as e:
            gridpp.apply_curve([0, 1], [[1, 2], [1, 2, 3]], gridpp.OneToOne, gridpp.OneToOne)

    def test_empty_fcst(self):
        """Check for empty result on empty input"""
        q = gridpp.apply_curve([], [[1, 2], [1, 2]], gridpp.OneToOne, gridpp.OneToOne)
        np.testing.assert_array_equal(q, [])

    def test_edge(self):
        """Check values on edge of curve"""
        x = [1, 2, 3]
        y = [2, 3, 4]
        curve = [x, y]
        policies = [gridpp.OneToOne, gridpp.Zero, gridpp.MeanSlope, gridpp.NearestSlope]
        for val in [1, 3]:
            for policy in policies:
                output = gridpp.apply_curve([val], curve, policy, policy)
                np.testing.assert_array_equal(output, [val+1])

    def test_extrapolation(self):
        """Check values outside curve"""
        x = [1, 2, 3]
        y = [2, 5, 6]
        curve = [x, y]
        policies = [gridpp.OneToOne, gridpp.Zero, gridpp.MeanSlope, gridpp.NearestSlope]
        np.testing.assert_array_equal(gridpp.apply_curve([0,4], curve, gridpp.OneToOne, gridpp.OneToOne), [1,7])
        np.testing.assert_array_equal(gridpp.apply_curve([0,4], curve, gridpp.Zero, gridpp.Zero), [2,6])
        np.testing.assert_array_equal(gridpp.apply_curve([0,4], curve, gridpp.MeanSlope, gridpp.MeanSlope), [0,8])
        np.testing.assert_array_equal(gridpp.apply_curve([0,4], curve, gridpp.NearestSlope, gridpp.NearestSlope), [-1,7])

if __name__ == '__main__':
    unittest.main()
