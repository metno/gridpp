import gridpp
import unittest
import numpy as np
import random

class Test(unittest.TestCase):
    def test_simple(self):
        base = np.expand_dims([0, 1, 2, 7, 15], 0)
        values = np.expand_dims([0, 1, 2, 1, 0], 0)
        halfwidth = 1
        min_num = 0
        min_range = 0
        default_gradient = -11
        gradient = gridpp.calc_gradient(base, values, gridpp.LinearRegression, halfwidth, min_num, min_range, default_gradient)
        np.testing.assert_array_almost_equal(gradient, [[1, 1, -0.064516, -0.151163, -1/8]])

    def test_small(self):
        """ Check when halfwidth is larger than the array """
        base = np.expand_dims([0, 1, 2], 0)
        values = np.expand_dims([0, 1, 2], 0)
        halfwidth = 5
        min_num = 0
        min_range = 0
        default_gradient = -11
        gradient = gridpp.calc_gradient(base, values, gridpp.LinearRegression, halfwidth, min_num, min_range, default_gradient)
        np.testing.assert_array_almost_equal(gradient, [[1, 1, 1]])

    def test_invalid_arguments(self):
        args = [gridpp.LinearRegression, 1, 0, 0, 0]
        with self.assertRaises(ValueError) as e:
            gridpp.calc_gradient(np.zeros([3, 2]), np.zeros([2, 3]), *args)

if __name__ == '__main__':
    unittest.main()
