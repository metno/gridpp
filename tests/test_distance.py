from __future__ import print_function
import unittest
import gridpp
import numpy as np


class Test(unittest.TestCase):
    def test_point_to_grid_geodetic(self):
        lons, lats = np.meshgrid([0, 1, 2], [0, 1])
        grid = gridpp.Grid(lats, lons)
        points = gridpp.Points([0, 0], [0, 0.6])
        # TODO
        # np.testing.assert_array_almost_equal(gridpp.distance(points, grid, 1), [[0, 44527.79], [111319.49, 119893.92], [222638.98,227046.33]])

    def test_point_to_grid_cartesian(self):
        lons, lats = np.meshgrid([0, 1000, 2000], [0, 1000])
        grid = gridpp.Grid(lats, lons, 0*lats, 0*lats, gridpp.Cartesian)
        points = gridpp.Points([0, 0], [0, 600], [0,0], [0,0], gridpp.Cartesian)
        np.testing.assert_array_almost_equal(gridpp.distance(points, grid, 1), [[0, 400, 1400], [1000, np.sqrt(1000**2 + 400**2), np.sqrt(1000**2 + 1400**2)]], 4)
        np.testing.assert_array_almost_equal(gridpp.distance(points, grid, 2), [[600, 1000, 2000], [np.sqrt(1000**2 + 600**2), np.sqrt(2) * 1000, np.sqrt(1000**2 + 2000**2)]], 4)
        np.testing.assert_array_almost_equal(gridpp.distance(points, grid, 10), [[600, 1000, 2000], [np.sqrt(1000**2 + 600**2), np.sqrt(2) * 1000, np.sqrt(1000**2 + 2000**2)]], 4)

    def test_grid_to_point_geodetic(self):
        lons, lats = np.meshgrid([0, 1, 2], [0, 1])
        grid = gridpp.Grid(lats, lons)
        points = gridpp.Points([0, 0], [0, 0.6])
        np.testing.assert_array_almost_equal(gridpp.distance(grid, points, 1), [0, 44528], 0)
        np.testing.assert_array_almost_equal(gridpp.distance(grid, points, 2), [111319.49, 66791.7], 0)
        np.testing.assert_array_almost_equal(gridpp.distance(grid, points, 10), [248907.83, 191514.84], 0)

    def test_grid_to_point_cartesian(self):
        lons, lats = np.meshgrid([0, 1000, 2000], [0, 1000])
        grid = gridpp.Grid(lats, lons, 0*lats, 0*lats, gridpp.Cartesian)
        points = gridpp.Points([0, 0], [0, 600], [0,0], [0,0], gridpp.Cartesian)
        np.testing.assert_array_almost_equal(gridpp.distance(grid, points, 1), [0, 400], 4)
        np.testing.assert_array_almost_equal(gridpp.distance(grid, points, 2), [1000, 600], 4)
        np.testing.assert_array_almost_equal(gridpp.distance(grid, points, 10), [np.sqrt(1000**2 + 2000**2), np.sqrt(1000**2 + 1400**2)], 4)


if __name__ == '__main__':
    unittest.main()
