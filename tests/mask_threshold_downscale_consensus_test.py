from __future__ import print_function
import unittest
import gridpp
import numpy as np


class Test(unittest.TestCase):
    def setUp(self):
        lons1, lats1 = np.meshgrid([10, 30], [50, 30])
        lons2, lats2 = np.meshgrid([5, 15, 25], [45, 35, 25])
        self.grid1 = gridpp.Grid(lats1, lons1)
        self.grid2 = gridpp.Grid(lats2, lons2)
        self.threshold_values = np.moveaxis(np.array([ [[-1., -1.],[-1., -1.]], 
                            [[0., 0.],  [0., 0.]], 
                            [[1., 1.],  [1., 1.]]
                            ]), 0, -1)
        self.thresholds = np.array([ [-2., -0.5, 0.5],
                                [0., 1., -1.], 
                                [2., 0.5, 0.]])
        self.valuestrue = np.moveaxis(np.array([ [[10., 5.],[3., 2.]], 
                            [[0., 1.],  [4., 0.]], 
                            [[3., 0.],  [0., 6.]]
                            ]), 0, -1)
        self.valuesfalse = np.moveaxis(np.array([ [[0., 0.],[0., 0.]], 
                            [[0., 0.],  [0., 0.]], 
                            [[0., 0.],  [0., 0.]]
                            ]), 0, -1)

    def test_mask_threshold_leq_downscale_consensus_mean(self):
        output = gridpp.mask_threshold_downscale_consensus(self.grid1, self.grid2, self.valuestrue, self.valuesfalse, self.threshold_values, self.thresholds, gridpp.Leq, gridpp.Mean)
        np.testing.assert_array_equal(output, np.array([[0., 3 + 1./3., 2.], [2 + 1./3., 2 + 1./3., 2./3.], [2 + 1./3., 2 + 1./3., 2./3.]], dtype=np.float32))

    def test_mask_threshold_leq_downscale_consensus_sum(self):
        output = gridpp.mask_threshold_downscale_consensus(self.grid1, self.grid2, self.valuestrue, self.valuesfalse, self.threshold_values, self.thresholds, gridpp.Leq, gridpp.Sum)
        np.testing.assert_array_equal(output, np.array([[0., 10., 6.], [7., 7., 2.], [7., 7., 2.]], dtype=np.float32))

    def test_mask_threshold_gt_downscale_consensus_median(self):
        output = gridpp.mask_threshold_downscale_consensus(self.grid1, self.grid2, self.valuestrue, self.valuesfalse, self.threshold_values, self.thresholds, gridpp.Gt, gridpp.Median)
        np.testing.assert_array_equal(output, np.array([[3., 0., 0.], [0., 0., 0.], [0., 0., 0.]], dtype=np.float32))

    def test_mask_threshold_lt_downscale_consensus_max(self):
        output = gridpp.mask_threshold_downscale_consensus(self.grid1, self.grid2, self.valuestrue, self.valuesfalse, self.threshold_values, self.thresholds, gridpp.Lt, gridpp.Max)
        np.testing.assert_array_equal(output, np.array([[0., 10., 5.], [3., 4., 0.], [4., 4., 2.]], dtype=np.float32))    

    def test_mask_threshold_geq_downscale_consensus_count(self):
        self.threshold_values[0, 1, 0] = np.NaN
        output = gridpp.mask_threshold_downscale_consensus(self.grid1, self.grid2, self.valuestrue, self.valuesfalse, self.threshold_values, self.thresholds, gridpp.Geq, gridpp.Count)
        np.testing.assert_array_equal(output, np.array([[3., 3., 2.], [3., 3., 3.], [3., 3., 3.]], dtype=np.float32))   

    