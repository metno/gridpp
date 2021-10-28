import gridpp
import unittest
import numpy as np 
import random

class Test(unittest.TestCase):
    def setUp(self):
        self.input_base = [[0, 0.5],[0.7,1.0]]
        self.input_values = [[10,15],[17,20]]

        self.input_base2 = [[0.9, 0.9],[1.0, 1.0]]
        self.input_values2 = [[12,14],[13,18]]

        self.input_base_nan = [[np.nan, 0.5],[0.7,1.0]]
        self.input_values_nan = [[10,15],[np.nan, 20]]
    """
    def test_size_of_inputs(self):
        np.testing.assert_array_equal(len(self.input_base),
            len(gridpp.neighbourhood_search(self.input_values, self.input_base, 1, 0, 0.95, 0.7, 1.0, 0.1)))

    def test_results(self):
        np.testing.assert_array_equal([[18.5, 18.5] , [18.5, 20]],
            gridpp.neighbourhood_search(self.input_values, self.input_base,1,0,0.95,0.7,1.0,0.1))

    def test_no_function_results(self):
        np.testing.assert_array_equal(self.input_values2,
            gridpp.neighbourhood_search(self.input_values2, self.input_base2,1,0,0.85,0.7,1.0,0.1))

    def test_nan_input(self):
        np.testing.assert_array_equal([[10, 20],[20, 20]], gridpp.neighbourhood_search(self.input_values_nan, self.input_base_nan ,1,0,0.95,0.7, 1.0, 0.1))

    """
if __name__ == '__main__':
    unittest.main()
