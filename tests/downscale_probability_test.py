from __future__ import print_function
import unittest
import gridpp
import numpy as np


class Test(unittest.TestCase):
    def setUp(self) -> None:
        lons1, lats1 = np.meshgrid([10, 30], [50, 30])
        lons2, lats2 = np.meshgrid([5, 15, 25], [45, 35, 25])
        self.grid1 = gridpp.Grid(lats1, lons1)
        self.grid2 = gridpp.Grid(lats2, lons2)
        self.values = np.array([ [[-1., -1.],[-1., -1.]], 
                            [[0., 0.],  [0., 0.]], 
                            [[1., 1.],  [1., 1.]]
                            ])
        self.thresholds = np.array([ [-2., -0.5, 0.5],
                                [0., 1., -1.], 
                                [2., 0.5, 0.]])

    def test_downscale_probability_leq(self):
        output = gridpp.downscale_probability(self.grid1, self.grid2, self.values, self.thresholds, gridpp.leq)
        np.testing.assert_array_equal(output, np.array([[0., 1./3., 2./3.], [2./3., 1., 1./3.], [1., 2./3., 2./3.]], dtype=np.float32))

    def test_downscale_probability_gt(self):
        output = gridpp.downscale_probability(self.grid1, self.grid2, self.values, self.thresholds, gridpp.gt)
        np.testing.assert_array_equal(output, np.array([[1., 2./3., 1./3.], [1./3., 0., 2./3.], [0., 1./3., 1./3.]], dtype=np.float32))

    def test_downscale_probability_geq(self):
        self.values[0, 1, 1] = np.NaN
        output = gridpp.downscale_probability(self.grid1, self.grid2, self.values, self.thresholds, gridpp.geq)
        np.testing.assert_array_equal(output, np.array([[1., 2./3., 1./3.], [2./3., 1./3., 1.], [0., 1./3., 1.]], dtype=np.float32))

    def test_downscale_probability_lt(self):
        self.values[:, 0, 0] = np.NaN
        output = gridpp.downscale_probability(self.grid1, self.grid2, self.values, self.thresholds, gridpp.lt)
        np.testing.assert_array_equal(output, np.array([[np.NaN, np.NaN, 2./3.], [1./3., 2./3., 0.], [1., 2./3., 1./3.]], dtype=np.float32))
