from __future__ import print_function
import unittest
import gridpp
import numpy as np


class Test(unittest.TestCase):
    def test_get_statistic(self):
        self.assertEqual(gridpp.get_statistic("mean"), gridpp.Mean)
        self.assertEqual(gridpp.get_statistic("min"), gridpp.Min)
        self.assertEqual(gridpp.get_statistic("max"), gridpp.Max)
        self.assertEqual(gridpp.get_statistic("count"), gridpp.Count)
        self.assertEqual(gridpp.get_statistic("median"), gridpp.Median)
        self.assertEqual(gridpp.get_statistic("quantile"), gridpp.Quantile)
        self.assertEqual(gridpp.get_statistic("std"), gridpp.Std)
        self.assertEqual(gridpp.get_statistic("sum"), gridpp.Sum)
        self.assertEqual(gridpp.get_statistic("randomchoice"), gridpp.RandomChoice)

    def test_unknown_statistic(self):
        self.assertEqual(gridpp.get_statistic("mean1"), gridpp.Unknown)

    """ Check that it doesn't cause any errors """
    def test_version(self):
        gridpp.version()

    def test_clock(self):
        time = gridpp.clock()
        self.assertTrue(time > 0)

    def test_is_valid(self):
        self.assertTrue(gridpp.is_valid(1))
        self.assertTrue(gridpp.is_valid(-1))
        self.assertTrue(gridpp.is_valid(-999))  # Check that the old missing value indicator is valid now
        self.assertFalse(gridpp.is_valid(np.nan))

    def test_calc_statistic_mean(self):
        self.assertEqual(gridpp.calc_statistic([0, 1, 2], gridpp.Mean), 1)
        self.assertEqual(gridpp.calc_statistic([0, 1, np.nan], gridpp.Mean), 0.5)
        self.assertEqual(gridpp.calc_statistic([np.nan, 1, np.nan], gridpp.Mean), 1)
        self.assertTrue(np.isnan(gridpp.calc_statistic([np.nan, np.nan, np.nan], gridpp.Mean)))
        self.assertTrue(np.isnan(gridpp.calc_statistic([], gridpp.Mean)))

    def test_calc_statistic_count(self):
        self.assertEqual(gridpp.calc_statistic([0, 1, 2], gridpp.Count), 3)
        self.assertEqual(gridpp.calc_statistic([0, 1, np.nan], gridpp.Count), 2)
        self.assertEqual(gridpp.calc_statistic([np.nan, 1, np.nan], gridpp.Count), 1)
        self.assertEqual(gridpp.calc_statistic([np.nan, np.nan, np.nan], gridpp.Count), 0)
        self.assertEqual(gridpp.calc_statistic([], gridpp.Count), 0)

    def test_calc_statistic_sum(self):
        self.assertEqual(gridpp.calc_statistic([0, 1, 2], gridpp.Sum), 3)
        self.assertEqual(gridpp.calc_statistic([0, 1, np.nan], gridpp.Sum), 1)
        self.assertEqual(gridpp.calc_statistic([np.nan, 1, np.nan], gridpp.Sum), 1)
        self.assertTrue(np.isnan(gridpp.calc_statistic([np.nan, np.nan, np.nan], gridpp.Sum)))
        self.assertTrue(np.isnan(gridpp.calc_statistic([], gridpp.Sum)))

    def test_calc_quantile(self):
        self.assertTrue(np.isnan(gridpp.calc_quantile([], 0)))
        self.assertEqual(gridpp.calc_quantile([0, 1, 2], 0), 0)
        self.assertEqual(gridpp.calc_quantile([0, 1, 2], 0.5), 1)
        self.assertEqual(gridpp.calc_quantile([0, 1, 2], 1), 2)
        self.assertEqual(gridpp.calc_quantile([0, np.nan, 2], 1), 2)
        self.assertEqual(gridpp.calc_quantile([0, np.nan, 2], 0), 0)
        self.assertEqual(gridpp.calc_quantile([0, np.nan, 2], 0.5), 1)
        for quantile in [0, 0.5, 1]:
            self.assertTrue(np.isnan(gridpp.calc_quantile([np.nan, np.nan, np.nan], quantile)))
            self.assertTrue(np.isnan(gridpp.calc_quantile([np.nan], quantile)))
        # BUG: This should work:
        # self.assertTrue(np.isnan(gridpp.calc_quantile([], 0.5)))

        self.assertEqual(gridpp.calc_quantile([[0, 1, 2]], 0), [0])
        self.assertEqual(gridpp.calc_quantile([[0, 1, 2]], 0.5), [1])
        self.assertEqual(gridpp.calc_quantile([[0, 1, 2]], 1), [2])
        self.assertEqual(gridpp.calc_quantile([[0, np.nan, 2]], 1), [2])
        self.assertEqual(gridpp.calc_quantile([[0, np.nan, 2]], 0), [0])
        self.assertEqual(gridpp.calc_quantile([[0, np.nan, 2]], 0.5), [1])
        quantile_of_nan_list = gridpp.calc_quantile([[np.nan, np.nan, np.nan]], 0.5)
        self.assertEqual(len(quantile_of_nan_list), 1)
        self.assertTrue(np.isnan(quantile_of_nan_list[0]))

    def test_calc_quantile_spatially_varying(self):
        Q = 5
        Y = 3
        X = 2
        input = np.reshape(np.arange(Y * X * Q), [Y, X, Q])
        levels = np.reshape([0.25, 0.6]*int(Y*X/2), [Y, X])
        output = gridpp.calc_quantile(input, levels)
        expected = np.reshape([1, 7.4, 11, 17.4, 21, 27.4], [Y, X])
        np.testing.assert_array_almost_equal(output, expected)

    def test_calc_quantile_spatially_varying_invalid_arguments(self):
        with self.assertRaises(ValueError) as e:
            # Dimension mismatch
            gridpp.calc_quantile(np.zeros([2, 3, 4]), np.zeros([1,3]))

    def test_calc_quantile_spatially_varying_empty(self):
        output = gridpp.calc_quantile(np.zeros([2, 0, 3]), np.zeros([2, 0]))
        np.testing.assert_array_almost_equal(output.shape, [0, 0])

        output = gridpp.calc_quantile(np.zeros([0, 2, 3]), np.zeros([0, 2]))
        np.testing.assert_array_almost_equal(output.shape, [0, 0])

        output = gridpp.calc_quantile(np.zeros([2, 2, 0]), np.zeros([2, 2]))
        np.testing.assert_array_almost_equal(output, np.nan*np.zeros([2, 2]))

    def test_calc_quantile_invalid_argument(self):
        quantiles = [1.1, -0.1]
        for quantile in quantiles:
            with self.assertRaises(ValueError) as e:
                gridpp.calc_quantile([0, 1, 2], quantile)
            with self.assertRaises(ValueError) as e:
                gridpp.calc_quantile([[0, 1, 2]], quantile)
            with self.assertRaises(ValueError) as e:
                gridpp.calc_quantile(np.zeros([2,3,4]), quantile*np.ones([2,3]))
        self.assertTrue(np.isnan(gridpp.calc_quantile([0, 1, 2], np.nan)))

    def test_calc_statistic_randomchoice(self):
        # since this is random, just check that we don't get unreasonable results
        for i in range(10):
            self.assertGreaterEqual(gridpp.calc_statistic([0, 1, 2], gridpp.RandomChoice), 0)
            self.assertLessEqual(gridpp.calc_statistic([0, 1, 2], gridpp.RandomChoice), 2)

    def test_calc_statistic_randomchoice_missing(self):
        ar = [1, np.nan, 2, 3, np.nan, np.nan]
        for i in range(10):
            output = gridpp.calc_statistic(ar, gridpp.RandomChoice)
            self.assertTrue(output in [1, 2, 3])

    def test_calc_statistic_randomchoice_most_missing(self):
        ar = np.nan * np.zeros([1000])
        ar[100] = 1
        self.assertEqual(gridpp.calc_statistic(ar, gridpp.RandomChoice), 1)

    def test_calc_statistic_randomchoice_only_missing(self):
        ar = np.nan * np.zeros([1000])
        self.assertTrue(np.isnan(gridpp.calc_statistic(ar, gridpp.RandomChoice)))

    def test_num_missing_values(self):
        self.assertEqual(gridpp.num_missing_values([[0, np.nan, 1, np.nan]]), 2)
        self.assertEqual(gridpp.num_missing_values([[np.nan, np.nan]]), 2)
        self.assertEqual(gridpp.num_missing_values([[0, 0, 1, 1]]), 0)
        self.assertEqual(gridpp.num_missing_values([[0, np.nan], [1, np.nan]]), 2)
        self.assertEqual(gridpp.num_missing_values([[np.nan, np.nan], [np.nan, np.nan]]), 4)
        self.assertEqual(gridpp.num_missing_values([[0, 0], [1, 1]]), 0)
        self.assertEqual(gridpp.num_missing_values([[]]), 0)

    def test_calc_statistics_2d(self):
        values = np.reshape(np.arange(9), [3, 3])
        output = gridpp.calc_statistic(values, gridpp.Mean)
        np.testing.assert_array_almost_equal(output, [1, 4, 7])

    def test_warning(self):
        gridpp.warning("test")

    def test_error(self):
        with self.assertRaises(RuntimeError) as e:
            gridpp.error("test")

    def test_get_index(self):
        self.assertEqual(2, gridpp.get_lower_index(1, [0, 0, 1, 1]))
        self.assertEqual(3, gridpp.get_upper_index(1, [0, 0, 1, 1]))
        self.assertEqual(0, gridpp.get_lower_index(0, [0, 0, 1, 1]))
        self.assertEqual(1, gridpp.get_upper_index(0, [0, 0, 1, 1]))
    def test_compatible_size_grid_vec2(self):
        lons, lats = np.meshgrid([0, 10, 20, 30], [30, 40, 50])
        grid = gridpp.Grid(lats, lons)  # 3 x 4
        self.assertFalse(gridpp.compatible_size(grid, np.zeros([2, 4])))
        self.assertTrue(gridpp.compatible_size(grid, np.zeros([3, 4])))

    def test_compatible_size_grid_vec3(self):
        lons, lats = np.meshgrid([0, 10, 20, 30], [30, 40, 50])
        grid = gridpp.Grid(lats, lons)  # 3 x 4
        self.assertFalse(gridpp.compatible_size(grid, np.zeros([3, 2, 4])))
        self.assertTrue(gridpp.compatible_size(grid, np.zeros([3, 3, 4])))

    def test_compatible_size_points_vec(self):
        points = gridpp.Points([0, 1, 2], [0, 1, 2])
        self.assertTrue(gridpp.compatible_size(points, [0, 0, 0]))

    def test_compatible_size_points_vec2(self):
        points = gridpp.Points([0, 1, 2], [0, 1, 2])
        self.assertTrue(gridpp.compatible_size(points, np.zeros([1, 3])))
        self.assertTrue(gridpp.compatible_size(points, np.zeros([2, 3])))

    def test_set_omp_threads(self):
        num = gridpp.get_omp_threads()
        if num == 0:
            gridpp.set_omp_threads(3)
            self.assertEqual(gridpp.get_omp_threads(), 0)
        else:
            gridpp.set_omp_threads(3)
            self.assertEqual(gridpp.get_omp_threads(), 3)

    def test_is_valid_lat(self):
        valid = [0, -90, 90]
        invalid = [-91, 91, np.nan]
        for v in valid:
            with self.subTest(v=v):
                self.assertTrue(gridpp.is_valid_lat(v, gridpp.Geodetic))
        for v in invalid:
            with self.subTest(v=v):
                self.assertFalse(gridpp.is_valid_lat(v, gridpp.Geodetic))

    def test_is_valid_lon(self):
        valid = [-1000, -180, 0, 180, 1000]
        invalid = [np.nan]
        for v in valid:
            with self.subTest(v=v):
                self.assertTrue(gridpp.is_valid_lon(v, gridpp.Geodetic))
        for v in invalid:
            with self.subTest(v=v):
                self.assertFalse(gridpp.is_valid_lon(v, gridpp.Geodetic))

    def test_set_debug_level(self):
        for level in [0, 1, 10]:
            gridpp.set_debug_level(level)
            self.assertEqual(level, gridpp.get_debug_level())

    def test_debug(self):
        gridpp.debug("test")
        gridpp.debug("")

    def test_future_deprecation_warning(self):
        gridpp.future_deprecation_warning("test")
        gridpp.future_deprecation_warning("")
        gridpp.future_deprecation_warning("test", "other")
        gridpp.future_deprecation_warning("", "other")
        gridpp.future_deprecation_warning("", "")

    def test_calc_even_quantiles(self):
        output = gridpp.calc_even_quantiles([1,2,3], 2)
        np.testing.assert_array_almost_equal(output, [1,3])

        # First and last
        output = gridpp.calc_even_quantiles(range(10), 2)
        np.testing.assert_array_almost_equal(output, [0,9])

        # Repeated first number
        output = gridpp.calc_even_quantiles([1,1,1,1,1,5, 10], 3)
        np.testing.assert_array_almost_equal(output, [1, 5, 10])

        output = gridpp.calc_even_quantiles([1,1,1,1,1,5, 10], 2)
        np.testing.assert_array_almost_equal(output, [1, 10])

        output = gridpp.calc_even_quantiles([1,1,1,1,4,5, 10], 3)
        np.testing.assert_array_almost_equal(output, [1, 4, 10])

        # Repeated numbers
        for num in [1, 2, 3]:
            with self.subTest(num=num):
                output = gridpp.calc_even_quantiles([1,1,1], num)
                np.testing.assert_array_almost_equal(output, [1])

        # Too little data
        output = gridpp.calc_even_quantiles([1,2,3], 3)
        np.testing.assert_array_almost_equal(output, [1,2,3])

        # Too little data with repeated
        output = gridpp.calc_even_quantiles([1,1,3], 3)
        np.testing.assert_array_almost_equal(output, [1,3])

        output = gridpp.calc_even_quantiles([1], 2)
        np.testing.assert_array_almost_equal(output, [1])

        # Empty arrays
        output = gridpp.calc_even_quantiles([1,2,3], 0)
        np.testing.assert_array_almost_equal(output, [])

        output = gridpp.calc_even_quantiles([], 0)
        np.testing.assert_array_almost_equal(output, [])

        output = gridpp.calc_even_quantiles([], 2)
        np.testing.assert_array_almost_equal(output, [])

if __name__ == '__main__':
    unittest.main()
