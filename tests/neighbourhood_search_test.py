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

        self.input_base_expand = [[0,0.5,0],[0,0.5,0],[0.75,1.0,0.75]]
        self.input_values_expand = [[10,15,10],[10,15,10],[17.5,20,17.5]]



    def get_apply_array(self, base, min_value, max_value):
        return (np.array(base) >= min_value) & (np.array(base) <= max_value)

    def test_size_of_inputs(self):
        apply_array = self.get_apply_array(self.input_base, 0, 0.95)
        np.testing.assert_array_equal(len(self.input_base),
            len(gridpp.neighbourhood_search(self.input_values, self.input_base, 1, 0.7, 1.0, 0.1, apply_array)))

    def test_results(self):
        apply_array = self.get_apply_array(self.input_base, 0, 0.95)
        np.testing.assert_array_equal([[18.5, 18.5] , [18.5, 20]],
            gridpp.neighbourhood_search(self.input_values, self.input_base, 1,0.7,1.0,0.1, apply_array))

    def test_no_apply_array(self):
        np.testing.assert_array_equal([[18.5, 18.5] , [18.5, 18.5]],
            gridpp.neighbourhood_search(self.input_values, self.input_base, 1,0.7,1.0,0.1))

    def test_no_function_results(self):
        apply_array = self.get_apply_array(self.input_base2, 0, 0.85)
        np.testing.assert_array_equal(self.input_values2,
            gridpp.neighbourhood_search(self.input_values2, self.input_base2,1,0.7,1.0,0.1, apply_array))

    def test_nan_input(self):
        apply_array = self.get_apply_array(self.input_base_nan, 0, 0.95)
        np.testing.assert_array_equal([[10, 20],[20, 20]],
                gridpp.neighbourhood_search(self.input_values_nan, self.input_base_nan, 1,0.7, 1.0, 0.1, apply_array))

    def test_simple(self):
        output = gridpp.neighbourhood_search([[0, 1, 2]], [[0.5, 0.5, 1]], 1, 0.7, 1, 0.1)
        np.testing.assert_array_equal(output, [[0, 2, 2]])

    def test_expand(self):
        apply_array = self.get_apply_array(self.input_base_expand, 0, 0.95)
        output = gridpp.neighbourhood_search(self.input_values_expand, self.input_base_expand, 1, 0.7, 1.0, 0.1, apply_array, True, 2.0/3, 1.0/2)
        np.testing.assert_array_almost_equal(output, [[15, 18 + 1.0/3, 15],[18.75, 18 + 1.0/3, 18.75],[18.75, 20, 18.75]])

    def test_expand_nan(self):
        input_base_expand_nan = self.input_base_expand
        input_base_expand_nan[1][1] = np.nan
        apply_array = self.get_apply_array(self.input_base_expand, 0, 0.95)
        output = gridpp.neighbourhood_search(self.input_values_expand, input_base_expand_nan, 1, 0.7, 1.0, 0.1, apply_array, True, 2.0/3, 1.0/2 ) 
        np.testing.assert_array_almost_equal(output, [[10,18 + 1.0/3,10],[10,15,10],[18.75,20,18.75]])

if __name__ == '__main__':
    unittest.main()


