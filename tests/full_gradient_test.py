import gridpp
import unittest
import numpy as np 
import random

class Test(unittest.TestCase):
    def setUp(self):
        gridpp.set_omp_threads(1)

        self.ilons, self.ilats = np.meshgrid([10,20], [40,50])
        self.olons, self.olats = np.meshgrid([10,15,20],[40,45,50])

        self.ilafs = np.ones([2,2])*0.5
        self.olafs = np.array([[0,1,0],[0,1,0],[0,1,0]])
        self.ielevs = np.zeros([2,2])
        self.oelevs = np.array([[0,100,0,],[0,200,0],[0,100,0]])
        self.ivalues = np.ones([2,2])*15
        self.elev_gradient = - np.ones(self.ilats.shape) / 100.0
        self.laf_gradient = np.ones(self.ilats.shape) * 10

        self.igrid = gridpp.Grid(self.ilats, self.ilons, self.ielevs, self.ilafs)
        self.ogrid = gridpp.Grid(self.olats, self.olons, self.oelevs, self.olafs)

    def test_grid_to_grid_laf_2d(self):
        self.elev_gradient[:,:] = 0

        output = gridpp.full_gradient(self.igrid, self.ogrid, self.ivalues, self.elev_gradient, self.laf_gradient)

        np.testing.assert_array_equal(output, [[10,20,10]]*3)

    def test_grid_to_grid_elev_2d(self):
        self.laf_gradient[:,:] = 0

        output = gridpp.full_gradient(self.igrid, self.ogrid, self.ivalues, self.elev_gradient, self.laf_gradient)

        np.testing.assert_array_equal(output, [[15,14,15], [15,13,15], [15,14,15]])

    def test_grid_to_grid_all_2d(self):
        output = gridpp.full_gradient(self.igrid, self.ogrid, self.ivalues, self.elev_gradient, self.laf_gradient)        

        np.testing.assert_array_equal(output, [[10,19,10], [10,18,10], [10,19,10]])

if __name__ == '__main__':
    unittest.main()