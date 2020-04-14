from __future__ import print_function
import unittest
import gridpp
import numpy as np


lats = [60, 60, 60, 60, 60, 70]
lons = [10,10.1,10.2,10.3,10.4, 10]


class NeighbourhoodTest(unittest.TestCase):
    def test_neighbourhood(self):
        """Simple check
        50  6  7  8
        40  3  4  5
        30  0  1  2
            0 10 20
        """
        lons, lats = np.meshgrid([0, 10, 20], [30, 40, 50])
        grid = gridpp.Grid(lats, lons)
        values = np.zeros(lons.shape)
        values[:] = np.reshape(range(9), lons.shape)
        points = gridpp.Points([25, 40, 55, 25, 40, 55, 25, 40, 55, 25, 40, 55, 25, 40, 55], [-1, -1, -1, 0, 0, 0, 10, 10, 10, 20, 20, 20, 21, 21, 21])
        output = gridpp.nearest(grid, points, values)
        self.assertEqual(output, (0, 3, 6, 0, 3, 6, 1, 4, 7, 2, 5, 8, 2, 5, 8))
        output = gridpp.bilinear(grid, points, values)
        self.assertEqual(output, (0, 3, 6, 0, 3, 6, 1, 4, 7, 2, 5, 8, 2, 5, 8))


if __name__ == '__main__':
    unittest.main()
