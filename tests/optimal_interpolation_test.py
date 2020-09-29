from __future__ import print_function
import unittest
import gridpp
import numpy as np


class Test(unittest.TestCase):
    """ Check if OI for a single observation works Put an obs at 2500m and background at 0, 2500,
    10000 m
    """
    def test_simple_1d(self):
        N = 3
        y = [[0, 0, 0]]
        x = [[0, 2500, 10000]]
        grid = gridpp.Grid(y, x, y, y, gridpp.Cartesian)
        points = gridpp.Points([0], [2500], [0], [0], gridpp.Cartesian)
        pratios = [0.1]
        structure = gridpp.BarnesStructure(2500)
        pobs = [1]
        background = np.zeros([1, N])
        pbackground = [0]
        max_points = 10
        output = gridpp.optimal_interpolation(grid, background, points, pobs, pratios, pbackground, structure, max_points)
        np.testing.assert_array_almost_equal(output, np.array([[np.exp(-0.5)/1.1, 1/1.1, np.exp(-0.5*9)/1.1]]))


if __name__ == '__main__':
    unittest.main()
