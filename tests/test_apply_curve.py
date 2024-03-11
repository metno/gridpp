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
            gridpp.apply_curve([0, 1], [1, 2, 3], [1, 2], gridpp.OneToOne, gridpp.OneToOne)
        with self.assertRaises(Exception) as e:
            gridpp.apply_curve([0, 1], [1, 2], [1, 2, 3], gridpp.OneToOne, gridpp.OneToOne)

    def test_empty_fcst(self):
        """Check for empty result on empty input"""
        q = gridpp.apply_curve([], [1, 2], [1, 2], gridpp.OneToOne, gridpp.OneToOne)
        np.testing.assert_array_equal(q, [])

    def test_edge(self):
        """Check values on edge of curve"""
        x = [1, 2, 3]
        y = [2, 3, 4]
        policies = [gridpp.OneToOne, gridpp.Zero, gridpp.MeanSlope, gridpp.NearestSlope]
        for val in [1, 3]:
            for policy in policies:
                output = gridpp.apply_curve([val], y, x, policy, policy)
                np.testing.assert_array_equal(output, [val+1])

    def test_extrapolation(self):
        """Check values outside curve"""
        x = [1, 2, 3]
        y = [2, 5, 6]
        policies = [gridpp.OneToOne, gridpp.Zero, gridpp.MeanSlope, gridpp.NearestSlope]
        np.testing.assert_array_equal(gridpp.apply_curve([0,4], y, x, gridpp.OneToOne, gridpp.OneToOne), [1,7])
        np.testing.assert_array_equal(gridpp.apply_curve([0,4], y, x, gridpp.Zero, gridpp.Zero), [2,6])
        np.testing.assert_array_equal(gridpp.apply_curve([0,4], y, x, gridpp.MeanSlope, gridpp.MeanSlope), [0,8])
        np.testing.assert_array_equal(gridpp.apply_curve([0,4], y, x, gridpp.NearestSlope, gridpp.NearestSlope), [-1,7])
        np.testing.assert_array_equal(gridpp.apply_curve([0,4], y, x, gridpp.Unchanged, gridpp.Unchanged), [0,4])

    def test_3d(self):
        curve_fcst = np.random.rand(3, 2, 4)
        curve_ref = np.random.rand(3, 2, 4)
        field = np.random.rand(3, 2)
        gridpp.apply_curve(field, curve_ref, curve_fcst, gridpp.OneToOne, gridpp.OneToOne)

    def test_all_extrapolation_policies(self):
        """Check that all policies work"""
        curve_fcst = np.random.rand(3, 2, 4)
        curve_ref = np.random.rand(3, 2, 4)
        field = np.random.rand(3, 2)
        for policy in [gridpp.OneToOne, gridpp.Zero, gridpp.NearestSlope, gridpp.MeanSlope, gridpp.Unchanged]:
            gridpp.apply_curve(field, curve_ref, curve_fcst, policy, policy)

    def test_invalid_extrapolation_policy(self):
        policy = -1
        curve_fcst = np.random.rand(3, 2, 4)
        curve_ref = np.random.rand(3, 2, 4)
        field = np.random.rand(3, 2)
        with self.assertRaises(ValueError) as e:
            gridpp.apply_curve(field, curve_ref, curve_fcst, policy, policy)


if __name__ == '__main__':
    unittest.main()
