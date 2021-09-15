from __future__ import print_function
import unittest
import gridpp
import numpy as np
import collections


class Test(unittest.TestCase):
    def test_no_obs(self):
        grid = gridpp.Points([0], [0])
        points = gridpp.Points([], [])
        pratios = []
        structure = gridpp.BarnesStructure(500)
        pobs = []
        background = np.zeros(grid.size())
        pbackground = []
        max_points = 10
        output0 = gridpp.optimal_interpolation(grid, background, points, pobs, pratios, pbackground, structure, max_points)
        np.testing.assert_almost_equal(output0, background)


if __name__ == '__main__':
    unittest.main()
