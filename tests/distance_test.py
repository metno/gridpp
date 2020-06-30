from __future__ import print_function
import unittest
import gridpp
import numpy as np


lats = [60, 60, 60, 60, 60, 70]
lons = [10,10.1,10.2,10.3,10.4, 10]


class DistanceTest(unittest.TestCase):
    def test_neighbourhood(self):
        """Simple check
         2  X
         1  X
         0  X
            0  1
        """
        for num_treads in [1, 2]:
            gridpp.set_omp_threads(num_treads)
            lons, lats = np.meshgrid([0, 1], [0, 1, 2])
            grid = gridpp.Grid(lats, lons)
            points = gridpp.Points([0, 1, 2], [0, 0, 0])
            radius = 120000
            output = gridpp.count(points, grid, radius)
            np.testing.assert_array_equal(output, ((2,1), (3,1), (2,1)))

    def test_one_row(self):
        gridpp.set_omp_threads(1)
        lons, lats = np.meshgrid([0], [30, 40, 50])
        grid = gridpp.Grid(lats, lons)
        values = np.zeros(lons.shape)
        values[:] = np.reshape(range(3), lons.shape)
        points = gridpp.Points([25, 40, 55, 25, 40, 55, 25, 40, 55, 25, 40, 55, 25, 40, 55], [-1, -1, -1, 0, 0, 0, 10, 10, 10, 20, 20, 20, 21, 21, 21])
        output = gridpp.nearest(grid, points, values)
        np.testing.assert_array_equal(output, (0, 1, 2, 0, 1, 2, 0, 1, 2, 0, 1, 2, 0, 1, 2))


if __name__ == '__main__':
    unittest.main()
