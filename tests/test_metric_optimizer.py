from __future__ import print_function
import unittest
import gridpp
import numpy as np
import time


class Test(unittest.TestCase):
    """ Test that the curve reproduces the reference values """
    def test_get_optimal_threshold(self):
        obs = [0, 1, 1.4, 1.6, 2]
        fcst = [1, 2, 2.4, 2.6, 3]
        threshold = 1.5
        metric = gridpp.Bias
        output = gridpp.get_optimal_threshold(obs, fcst, threshold, metric)
        self.assertTrue(output > 2.4 and output < 2.6)

    def test_metric_optimizer_curve(self):
        # Check that a quadratic relationship gets fixed
        obs = np.linspace(0, 10, 101)
        fcst = obs ** 2
        thresholds = [1, 2, 3, 4]
        for metric in [gridpp.Pc, gridpp.Bias]:
            with self.subTest(metric=metric):
                x, y = gridpp.metric_optimizer_curve(obs, fcst, thresholds, metric)
                expected = y**2
                np.testing.assert_array_almost_equal(expected, x, 0)

    def test_calc_score(self):
        # Use verif package to check results
        obs = [1, 2, 3]
        fcst = [2, 1, 3]
        thresholds = [-1, 1.5, 2.5, 3]
        expected = {
                gridpp.Bias: [1, 1, 1, 1],
                gridpp.Pc: [1, 0.333, 1, 1],
                gridpp.Ets: [np.nan, -0.2, 1, np.nan],
                gridpp.Kss: [np.nan, -0.5, 1, np.nan],
                gridpp.Hss: [np.nan, -0.5, 1, np.nan],
                gridpp.Ts: [1, 0.333, 1, np.nan]
            }
        for metric in expected.keys():
            for t, threshold in enumerate(thresholds):
                with self.subTest(metric=metric, threshold=threshold):
                    output = gridpp.calc_score(obs, fcst, threshold, metric)
                    np.testing.assert_almost_equal(expected[metric][t], output, 3)


if __name__ == '__main__':
    unittest.main()
