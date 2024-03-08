from __future__ import print_function
import unittest
import gridpp
import numpy as np


class Test(unittest.TestCase):
    def test_attributes(self):
        """Test that lat, lons, etc are set"""
        grid = gridpp.Grid([[0, 1], [1, 2]], [[3, 4], [4, 5]], [[6, 7], [7, 8]], [[0.1, 0.2], [0.1, 0.2]])
        np.testing.assert_array_almost_equal(grid.get_lats(), [[0, 1], [1, 2]])
        np.testing.assert_array_almost_equal(grid.get_lons(), [[3, 4], [4, 5]])
        np.testing.assert_array_almost_equal(grid.get_elevs(), [[6, 7], [7, 8]])
        np.testing.assert_array_almost_equal(grid.get_lafs(), [[0.1, 0.2], [0.1, 0.2]])

    def test_get_neighbours_with_distance(self):
        grid = gridpp.Grid([[0, 0, 0], [1000, 1000, 1000]], [[0, 1000, 2000], [0, 1000, 2000]],
                np.zeros([3, 2]), np.zeros([3, 2]), gridpp.Cartesian)
        indices, distances = grid.get_neighbours_with_distance(0, 0, 1500)
        self.assertEqual(len(indices), 4)
        np.testing.assert_array_almost_equal(distances, [0, 1000, 1000, np.sqrt(2)*1000], 4)

    def test_get_box(self):
        # 2x3 grid
        lats = [[0, 0, 0], [1, 1, 1]]
        lons = [[0, 1, 2], [0.25, 1.25, 2.25]]
        grid = gridpp.Grid(lats, lons)
        # Check that the point lat=0.4,lon=1.25 is in between lat[0:1] and lon[1:2]
        np.testing.assert_array_equal(grid.get_box(0.4, 1.25), [True, 0, 1, 1, 2])
        np.testing.assert_array_equal(grid.get_box(0.4, 1), [True, 0, 0, 1, 1])

    def test_get_box_almost_empty(self):
        """Check case where grid is 1x2 and no box exists"""
        grid = gridpp.Grid([[0,0]], [[0, 1]])
        self.assertFalse(grid.get_box(0.4, 1.25)[0])

    def test_get_box_empty(self):
        """Check case where grid is empty"""
        grid = gridpp.Grid()
        self.assertFalse(grid.get_box(0.4, 1.25)[0])

    def test_size(self):
        np.testing.assert_array_equal([2, 3], gridpp.Grid([[0, 0, 0], [1, 1, 1]], [[0, 1, 2], [0, 1, 2]]).size())
        np.testing.assert_array_equal([0, 0], gridpp.Grid([[], []], [[], []]).size())

    def test_empty(self):
        grid = gridpp.Grid()
        np.testing.assert_array_equal(grid.get_nearest_neighbour(0, 0), [])
        np.testing.assert_array_equal(grid.get_neighbours(0, 0, 1000), np.zeros([0, 0]))
        self.assertEqual(grid.get_num_neighbours(0, 0, 1000), 0)
        np.testing.assert_array_equal(grid.get_lats(), np.zeros([0, 0]))
        np.testing.assert_array_equal(grid.get_lons(), np.zeros([0, 0]))
        np.testing.assert_array_equal(grid.get_elevs(), np.zeros([0, 0]))
        np.testing.assert_array_equal(grid.get_lafs(), np.zeros([0, 0]))
        np.testing.assert_array_equal(grid.size(), [0, 0])

    def test_empty2(self):
        grid = gridpp.Grid(np.zeros([0, 0]), np.zeros([0, 0]))
        np.testing.assert_array_equal(grid.get_nearest_neighbour(0, 0), [])
        np.testing.assert_array_equal(grid.get_neighbours(0, 0, 1000), np.zeros([0, 0]))
        self.assertEqual(grid.get_num_neighbours(0, 0, 1000), 0)
        np.testing.assert_array_equal(grid.get_lats(), np.zeros([0, 0]))
        np.testing.assert_array_equal(grid.get_lons(), np.zeros([0, 0]))
        np.testing.assert_array_equal(grid.get_elevs(), np.zeros([0, 0]))
        np.testing.assert_array_equal(grid.get_lafs(), np.zeros([0, 0]))
        np.testing.assert_array_equal(grid.size(), [0, 0])

    def test_get_point(self):
        """Check that the right values are retrieved"""
        # 2x3 grid
        lats = [[0, 0, 0], [1, 1, 1]]
        lons = [[0, 1, 2], [0.25, 1.25, 2.25]]
        elevs = [[0, 1, 2], [3, 4, 5]]
        lafs = [[0, 1, 2], [3, 4, 5]]

        grid = gridpp.Grid(lats, lons, elevs, lafs)

        # Retrieve values for the second row (of 2), first column (of 3)
        p = grid.get_point(1, 0)
        self.assertEqual(p.lat, 1)
        self.assertEqual(p.lon, 0.25)
        self.assertEqual(p.elev, 3)
        self.assertEqual(p.laf, 3)

        s, x, y, z = gridpp.convert_coordinates(1, 0.25, gridpp.Geodetic)
        self.assertEqual(p.x, x)
        self.assertEqual(p.y, y)
        self.assertEqual(p.z, z)

if __name__ == '__main__':
    unittest.main()
