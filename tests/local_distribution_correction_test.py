from __future__ import print_function
import unittest
import gridpp
import numpy as np


class Test(unittest.TestCase):
    def test_parallel(self):
        """Check that the parallel version gives the same results"""
        lons, lats = np.meshgrid(np.linspace(0, 100000, 100), np.linspace(0, 100000, 100))
        grid = gridpp.Grid(lats, lons, 0*lats, 0*lons, gridpp.Cartesian)

        lats = np.random.rand(100) * 100000
        lons = np.random.rand(100) * 100000
        points = gridpp.Points(lats, lons, 0*lats, 0*lons, gridpp.Cartesian)
        a = np.random.rand(100, 100)
        b = np.random.rand(1000)
        c = np.random.rand(1000)
        structure = gridpp.BarnesStructure(2500)
        values = dict()
        for ncores in [1, 8]:
            gridpp.set_omp_threads(ncores)
            values[ncores] = gridpp.local_distribution_correction(grid, a,
                points, b, c, structure, 0.1, 0.9, 5)

        np.testing.assert_array_equal(values[1], values[8])


if __name__ == '__main__':
    unittest.main()
