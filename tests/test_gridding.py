from __future__ import print_function
import unittest
import gridpp
import numpy as np


def get_shape(grid):
    shape = grid.size()
    if isinstance(shape, int):
        shape = [shape,]
    return shape

class Test(unittest.TestCase):
    def setUp(self):
        y, x = np.meshgrid(np.linspace(0, 1, 2), np.linspace(0, 1, 3))
        self.grid = gridpp.Grid(y, x, 0*y, 0*y, gridpp.Cartesian)
        self.grid_as_points = self.grid.to_points()
        self.points = gridpp.Points([-0.2, 0.5, 1], [-0.2, 0.5, 1], [0, 0, 0], [0, 0, 0], gridpp.Cartesian)
        self.values = [1, 2, 3]
        self.radius = 0.6
        self.statistic = gridpp.Sum
        self.min_num = 0

    def test_min_num(self):
        expected_list = {
                0: [[1, np.nan], [2, 5], [np.nan, 3]],
                1: [[1, np.nan], [2, 5], [np.nan, 3]],
                2: [[np.nan, np.nan], [np.nan, 5], [np.nan, np.nan]],
            }

        for grid in [self.grid, self.grid_as_points]:
            for min_num, expected in expected_list.items():
                with self.subTest(grid=grid, min_num=min_num):
                    output = gridpp.gridding(grid, self.points, self.values, self.radius, min_num, self.statistic)
                    np.testing.assert_array_almost_equal(output.shape, get_shape(grid))
                    np.testing.assert_array_almost_equal(output.flatten(), np.array(expected).flatten())

    def test_statistic(self):
        expected_list = {
                gridpp.Sum:[[1, np.nan], [2, 5], [np.nan, 3]],
                gridpp.Mean: [[1, np.nan], [2, 2.5], [np.nan, 3]],
                gridpp.Count: [[1, 0], [1, 2], [0, 1]],
            }

        for grid in [self.grid, self.grid_as_points]:
            for statistic, expected in expected_list.items():
                with self.subTest(grid=grid, statistic=statistic):
                    output = gridpp.gridding(grid, self.points, self.values, self.radius, self.min_num, statistic)
                    np.testing.assert_array_almost_equal(output.shape, get_shape(grid))
                    np.testing.assert_array_almost_equal(output.flatten(), np.array(expected).flatten())

    def test_radius(self):
        expected_list = {
                0.001:[[np.nan, np.nan], [np.nan, np.nan], [np.nan, 3]],
                0.6: [[1, np.nan], [2, 5], [np.nan, 3]],
                10: [[6, 6], [6, 6], [6, 6]],
            }

        for grid in [self.grid, self.grid_as_points]:
            for radius, expected in expected_list.items():
                with self.subTest(grid=grid, radius=radius):
                    output = gridpp.gridding(grid, self.points, self.values, radius, self.min_num, self.statistic)
                    np.testing.assert_array_almost_equal(output.shape, get_shape(grid))
                    np.testing.assert_array_almost_equal(output.flatten(), np.array(expected).flatten())

    def test_invalid_arguments(self):
        # Values don't match point size
        for grid in [self.grid, self.grid_as_points]:
            with self.subTest(grid=grid):
                with self.assertRaises(ValueError) as e:
                    output = gridpp.gridding(grid, self.points, [0], self.radius, self.min_num, self.statistic)

    def test_empty_input(self):
        # Empty points should give all nans
        points = gridpp.Points([], [], [], [], gridpp.Cartesian)
        values = []
        expected = np.nan * np.zeros(self.grid.size())
        for grid in [self.grid, self.grid_as_points]:
            for statistic in [gridpp.Sum, gridpp.Mean]:
                output = gridpp.gridding(grid, points, values, self.radius, self.min_num, statistic)
                np.testing.assert_array_almost_equal(output.shape, get_shape(grid))
                np.testing.assert_array_almost_equal(output.flatten(), np.array(expected).flatten())

        # Count should be 0
        output = gridpp.gridding(self.grid, points, values, self.radius, self.min_num, gridpp.Count)
        expected = np.zeros(self.grid.size())
        np.testing.assert_array_almost_equal(output, expected)

    def test_empty_grid(self):
        grid = gridpp.Grid([[]], [[]], [[]], [[]], gridpp.Cartesian)
        for statistic in [gridpp.Sum, gridpp.Mean, gridpp.Count]:
            output = gridpp.gridding(grid, self.points, self.values, self.radius, self.min_num, statistic)
            np.testing.assert_array_almost_equal(output, np.zeros([0,0]))

    def test_empty_grid_as_points(self):
        grid_as_points = gridpp.Points([], [], [], [], gridpp.Cartesian)
        for statistic in [gridpp.Sum, gridpp.Mean, gridpp.Count]:
            output = gridpp.gridding(grid_as_points, self.points, self.values, self.radius, self.min_num, statistic)
            np.testing.assert_array_almost_equal(output, np.zeros([0]))

    def test_invalid_radius(self):
        radii = [-1, np.nan]
        for grid in [self.grid, self.grid_as_points]:
            for radius in radii:
                with self.subTest(grid=grid, radius=radius):
                    with self.assertRaises(ValueError) as e:
                        output = gridpp.gridding(grid, self.points, self.values, radius, self.min_num, self.statistic)

    def test_invalid_min_num(self):
        min_num = -1  # Can't check np.nan, since this is not an int, and will fail with a TypeError
        for grid in [self.grid, self.grid_as_points]:
            with self.subTest(grid=grid):
                with self.assertRaises(ValueError) as e:
                    output = gridpp.gridding(grid, self.points, self.values, self.radius, min_num, self.statistic)


if __name__ == '__main__':
    unittest.main()
