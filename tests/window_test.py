import gridpp
import unittest
import numpy as np 
import random

class Test(unittest.TestCase):
    def setUp(self):
        self.inputs = np.ones([5,5], float)
        for i, value in enumerate(self.inputs[:,0]):
            for j, value in enumerate(self.inputs[0,:]):
                self.inputs[i,j] = i + j

        self.inputs_nan = self.inputs.copy()
        self.inputs_nan[3,3] = np.nan

        self.small_inputs = np.ones([2,2], float)

    def test_sum(self):
        output = gridpp.window(self.inputs, 3, gridpp.Sum, False, False, False)
        np.testing.assert_array_equal(output, [[1,3,6,9, 7],[3,6,9,12,9],[5,9,12,15,11],[7,12,15,18,13],[9,15,18,21,15]])

    def test_mean(self):
        output = gridpp.window(self.inputs, 3, gridpp.Mean, False, False, False)
        np.testing.assert_array_equal(output, [[0.5, 1, 2, 3, 3.5],[1.5, 2, 3, 4, 4.5],[2.5, 3, 4, 5, 5.5],[3.5, 4, 5, 6, 6.5],[4.5, 5,6,7,7.5]])

    def test_min(self):
        output = gridpp.window(self.inputs, 3, gridpp.Min, False, False, False)
        np.testing.assert_array_equal(output, [[0,0,1,2,3],[1,1,2,3,4],[2,2,3,4,5],[3,3,4,5,6],[4,4,5,6,7]])
    
    def test_max(self):
        output = gridpp.window(self.inputs, 3, gridpp.Max, False, False, False) 
        np.testing.assert_array_equal(output, [[1,2,3,4,4],[2,3,4,5,5],[3,4,5,6,6],[4,5,6,7,7],[5,6,7,8,8]])    

    def test_sum_before(self):
        output = gridpp.window(self.inputs, 3, gridpp.Sum, True, False, False)
        np.testing.assert_array_equal(output, [[0,1,3,6,9],[1,3,6,9,12],[2,5,9,12,15],[3,7,12,15,18],[4,9,15,18,21]])

    def test_sum_missing_edge(self):
        output = gridpp.window(self.inputs, 3, gridpp.Sum, True, False, True)
        np.testing.assert_array_equal(output, [[np.nan, np.nan,3,6,9],[np.nan, np.nan,6,9,12],[np.nan, np.nan,9,12,15],[np.nan, np.nan,12,15,18],[np.nan, np.nan,15,18,21]])
    
    def test_sum_keep_missing(self):
        output = gridpp.window(self.inputs_nan, 3, gridpp.Sum, True, True, False)
        np.testing.assert_array_equal(output, [[0,1,3,6,9],[1,3,6,9,12],[2,5,9,12,15],[3,7,12,np.nan,np.nan],[4,9,15,18,21]])

    def test_edge_case(self):
        output = gridpp.window(self.small_inputs, 5, gridpp.Sum, False, False, False)
        np.testing.assert_array_equal(output, [[2,2],[2,2]])

    def test_edge_case2(self):
        output = gridpp.window(self.small_inputs, 5, gridpp.Sum, False, False, True)
        np.testing.assert_array_equal(output, [[np.nan, np.nan],[np.nan, np.nan]])

if __name__ == '__main__':
    unittest.main()


#if window is beigger than array () # Edge cases