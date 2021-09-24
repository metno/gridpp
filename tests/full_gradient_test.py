import gridpp
import unittest
import numpy as np 
import random

class Test(unittest.TestCase):
    def setUp(self):
        return 

    def test_grid_to_grid_laf(self):
        """Tests that laf_gradient creates correct output independantly"""
        gridpp.set_omp_threads(1)
        #initialize
        ilons, ilats = np.meshgrid([10,20], [40,50])
        olons, olats = np.meshgrid([10,15,20],[40,45,50])
        ilafs = np.zeros([2,2])
        olafs = np.zeros([3,3])
        ielevs = np.zeros([2,2])
        oelevs = np.zeros([3,3])
        ivalues = np.zeros([2,2])
        ivalues[:,:] = 13.3
        ilafs[:,:] = 0.33
        olafs[:,:] = [[0,1,0],[0,1,0],[0,1,0]]

        igrid = gridpp.Grid(ilats, ilons, ielevs, ilafs)
        ogrid = gridpp.Grid(olats, olons, oelevs, olafs)

        elev_gradient = np.zeros(ilats.shape)
        laf_gradient = np.ones(ilats.shape) * 10

        output = gridpp.full_gradient(igrid, ogrid, ivalues, elev_gradient, laf_gradient)

        np.testing.assert_array_equal(output, [[10,20,10]]*3)

    def test_grid_to_grid_elev(self):
        """ Tests that elev_gradient creates correct output"""
        gridpp.set_omp_threads(1)
        #initialize
        ilons, ilats = np.meshgrid([10,20], [40,50])
        olons, olats = np.meshgrid([10,15,20],[40,45,50])
        ilafs = np.zeros([2,2])
        olafs = np.zeros([3,3])
        ielevs = np.zeros([2,2])
        oelevs = np.zeros([3,3])
        ivalues = np.ones([2,2]) * 10

        oelevs[:,:] = [[0,100,0,],[0,200,0],[0,100,0]]

        igrid = gridpp.Grid(ilats, ilons, ielevs, ilafs)
        ogrid = gridpp.Grid(olats, olons, oelevs, olafs)

        elev_gradient = - np.ones(ilats.shape) / 100.0
        laf_gradient = np.zeros(ilats.shape)

        output = gridpp.full_gradient(igrid, ogrid, ivalues, elev_gradient, laf_gradient)
        
        np.testing.assert_array_equal(output, [[10,9,10], [10,8,10], [10,9,10]])

    def test_grid_to_grid_all(self):
        """ Tests that elev_gradient and laf_gradient creates correct output"""
        gridpp.set_omp_threads(1)
        #initialize
        ilons, ilats = np.meshgrid([10,20], [40,50])
        olons, olats = np.meshgrid([10,15,20],[40,45,50])
        ilafs = np.zeros([2,2])
        olafs = np.zeros([3,3])
        ielevs = np.zeros([2,2])
        oelevs = np.zeros([3,3])
        ivalues = np.zeros([2,2])
        ivalues[:,:] = 15
        ilafs[:,:] = 0.5
        ielevs[:,:] = 0
        olafs[:,:] = [[0,1,0],[0,1,0],[0,1,0]]
        oelevs[:,:] = [[0,100,0,],[0,200,0],[0,100,0]]

        igrid = gridpp.Grid(ilats, ilons, ielevs, ilafs)
        ogrid = gridpp.Grid(olats, olons, oelevs, olafs)

        elev_gradient = - np.ones(ilats.shape) / 100.0
        laf_gradient = np.ones(ilats.shape) * 10

        output = gridpp.full_gradient(igrid, ogrid, ivalues, elev_gradient, laf_gradient)
        
        np.testing.assert_array_equal(output, [[10,19,10], [10,18,10], [10,19,10]])

if __name__ == '__main__':
    unittest.main()