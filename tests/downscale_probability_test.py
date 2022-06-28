from __future__ import print_function
import unittest
import gridpp
import numpy as np


class Test(unittest.TestCase):
    def test_downscale_probability_leq(self):
        lons1, lats1 = np.meshgrid([10, 30], [50, 30])
        lons2, lats2 = np.meshgrid([5, 15, 25], [45, 35, 25])
        grid1 = gridpp.Grid(lats1, lons1)
        grid2 = gridpp.Grid(lats2, lons2)
        values = np.array([ [[-1., -1.],[-1., -1.]], 
                            [[0., 0.],  [0., 0.]], 
                            [[1., 1.],  [1., 1.]]
                            ])
        thresholds = np.array([ [-2., -0.5, 0.5],
                                [0., 1., -1.], 
                                [2., 0.5, 0.]])
        output = gridpp.downscale_probability(grid1, grid2, values, thresholds, gridpp.leq)
        np.testing.assert_array_equal(output, np.array([[0., 1./3., 2./3.], [2./3., 1., 1./3.], [1., 2./3., 2./3.]], dtype=np.float32))
