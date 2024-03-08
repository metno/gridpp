from __future__ import print_function
import unittest
import gridpp
import numpy as np
import collections


class Test(unittest.TestCase):
    def test_simple(self):
        N = 11
        x = np.linspace(0, 10000, N)
        y = np.linspace(0, 10000, N)
        xx, yy = np.meshgrid(x, y)
        grid = gridpp.Grid(xx, yy, 0*xx, 0*xx, gridpp.Cartesian)
        points = gridpp.Points([3000, 5000], [10000, 5000], [0, 0], [0, 0], gridpp.Cartesian)
        obs = [-10, 10]
        background = np.zeros([len(x), len(x)])
        max_elev_diff = 200
        half_width = [1, 1]
        output = gridpp.doping_square(grid, background, points, obs, half_width, max_elev_diff)

        # Expect a field of 0s, with one square of 10s at 5,5 and one (half) square of -10 at 10,3
        expected = np.zeros([N, N])
        expected[4:7, 4:7] = 10
        expected[9:11, 2:5] = -10
        np.testing.assert_array_almost_equal(output, expected)


if __name__ == '__main__':
    unittest.main()
