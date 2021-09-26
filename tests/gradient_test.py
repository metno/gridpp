import gridpp
import unittest
import numpy as np 
import random

class Test(unittest.TestCase):	
    def setUp(self):
        self.input_base = [[0 ,0 , 0.2, 0.4], [0, 0.2, 0.4, 0.6], [0.2, 0.4,0.6, 0.8],[0.4, 0.6, 0.8, 1.0]]
        #self.input_values = [[11, 11, 12, 14], [11, 12, 14, 18],[12, 14, 18, 26], [14, 18, 26, 42]] 
        self.input_values = [[10,10,12,14], [10,12,14,16],[12,14,16,18],[14,16,18,20]] 

    def test_size_of_inputs(self):
        """ Checks that the input size of input and output gradient matrix are at the same size""" 
        np.testing.assert_array_equal(len(self.input_base), len( gridpp.calc_gradient(self.input_base, self.input_values, 3))) 

        #np.testing.assert_array_equal(self.input_base.shape, gridpp.calc_gradient(self.input_base, self.input_values, 3).shape)
        #np.testing.assert_array_equal(self.input_values.shape, gridpp.calc_gradient(self.input_base, self.input_values, 3).shape)

    def test_gradient(self):
        """ tests that inputs give the correct gradient"""
        self.input_values = [[10,10,12,14], [10,12,14,16],[12,14,16,18],[14,16,18,20]]
        output = gridpp.calc_gradient(self.input_base, self.input_values, 3)
        np.testing.assert_array_almost_equal(output, np.ones(output.shape)*10, decimal=3)


    def test_nan_input(self):
        """ Check gradient of inputs, where bases and values are nan, output should be equal to default gradient"""
        input_base =  [[np.nan, np.nan, np.nan, np.nan]]*4
        default_gradient = -1
        np.testing.assert_array_almost_equal(gridpp.calc_gradient(input_base,self.input_values,
            3, 0, 0, default_gradient), default_gradient * np.ones(np.array(self.input_values).shape))

    def test_non_changing_base(self):
        """ If an matrix contains a field with no changing base numberrs, test the output is default gradient"""
        self.input_base = [[0,0,0,0],[0,0,0,0],[0,0,0,0],[0,0,0,0]]
        default_gradient = -1
        output = gridpp.calc_gradient(self.input_base, self.input_values, 3, 0, 0, default_gradient)
        np.testing.assert_equal(output[0,0], default_gradient)

    def test_empty_input(self):
        """ Check empty input returns and empty output. 
            Requires """
        return

    def test_nan_base(self):
        return

    def test_nan_values(self):
        return

    def blanding_nan_values(self):
        return 

    
    def test_zero_halfwidth(self):
        with self.assertRaises(ValueError) as e:
            gridpp.calc_gradient(self.input_base, self.input_values, -1)
    
    
    def test_negative_halfwidth(self):
        return


if __name__ == '__main__':
	unittest.main()


"""
Notes

self.input_base = [[0, 0, 0], [1, 1, 1], [4, 4, 4]]
self.input_base[0:2, 0:2] = np.nan

input size is null, return empty vector

check test_invalid arguments in neighbourhood_test.py line 31

use test with squared temperatures 
"""
