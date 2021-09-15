from __future__ import print_function
import unittest
import gridpp
import numpy as np
import collections


class Test(unittest.TestCase):
    def test_no_obs(self):
        """ Check that we get the same field back if there are no observations """
        grid = gridpp.Points([0], [0])
        E = 3
        points = gridpp.Points([], [])
        psigmas = []
        structure = gridpp.BarnesStructure(500000)
        pobs = []
        background = np.zeros([grid.size(), E])
        pbackground = np.zeros([0, E])
        max_points = 10
        output0 = gridpp.optimal_interpolation_ensi(grid, background, points, pobs, psigmas, pbackground, structure, max_points)
        np.testing.assert_almost_equal(output0, background)

    def test_some_missing_obs(self):
        """ Check that if one observation is missing, that the whole output isn't nan """
        grid = gridpp.Points([0], [0])
        E = 3
        points = gridpp.Points([0, 0.1], [0, 0.1])
        psigmas = [1, 1]
        structure = gridpp.BarnesStructure(500000)
        pobs = [np.nan, 0]
        background = np.zeros([grid.size(), E])
        pbackground = np.zeros([points.size(), E])
        max_points = 10
        output0 = gridpp.optimal_interpolation_ensi(grid, background, points, pobs, psigmas, pbackground, structure, max_points)
        np.testing.assert_almost_equal(output0, background)


if __name__ == '__main__':
    unittest.main()
