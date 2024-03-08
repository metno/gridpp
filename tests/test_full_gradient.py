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

        self.ivalues3d = np.ones([5,2,2])*15
        self.elev_gradient3d = - np.ones([5,2,2]) / 100.0 
        self.laf_gradient3d = np.ones([5,2,2]) * 10

        self.igrid = gridpp.Grid(self.ilats, self.ilons, self.ielevs, self.ilafs)
        self.ogrid = gridpp.Grid(self.olats, self.olons, self.oelevs, self.olafs)

        self.ipoints = gridpp.Points([10, 10, 20, 20], [45, 50, 45, 50], [0,0,0,0], [0.5, 0.5, 0.5, 0.5])
        self.opoints = gridpp.Points([10,10,10, 15, 15, 15, 20, 20, 20], [40, 45, 50, 40, 45, 50, 40, 45, 50], [0,100,0,0,200,0,0,100,0], [0,1,0,0,1,0,0,1,0])

    def test_grid_to_grid_all_2d(self):
        output = gridpp.full_gradient(self.igrid, self.ogrid, self.ivalues, self.elev_gradient, self.laf_gradient)        
        np.testing.assert_array_equal(output, [[10,19,10], [10,18,10], [10,19,10]])

    def test_grid_to_point_all_1d(self):
        output = gridpp.full_gradient(self.igrid, self.opoints, self.ivalues, self.elev_gradient, self.laf_gradient) 
        np.testing.assert_array_equal(output, [10,19,10,10,18,10,10,19,10])

        """Testing output using LAF parameter (elevation excluded)"""
    def test_grid_to_grid_laf_2d(self):
        self.elev_gradient[:,:] = 0
        output = gridpp.full_gradient(self.igrid, self.ogrid, self.ivalues, self.elev_gradient, self.laf_gradient)
        np.testing.assert_array_equal(output, [[10,20,10]]*3)

    def test_grid_to_grid_elev_2d(self):
        self.laf_gradient[:,:] = 0
        output = gridpp.full_gradient(self.igrid, self.ogrid, self.ivalues, self.elev_gradient, self.laf_gradient)
        np.testing.assert_array_equal(output, [[15,14,15], [15,13,15], [15,14,15]])

        """ Testing output using Elevation paramter (LAF excluded)"""
    def test_grid_to_point_laf_1d(self):       
        self.elev_gradient[:,:] = 0
        output = gridpp.full_gradient(self.igrid, self.opoints, self.ivalues, self.elev_gradient, self.laf_gradient) 
        np.testing.assert_array_equal(output, [10,20,10,10,20,10,10,20,10])

    def test_grid_to_point_elev_1d(self):
        self.laf_gradient[:,:] = 0 
        output = gridpp.full_gradient(self.igrid, self.opoints, self.ivalues, self.elev_gradient, self.laf_gradient) 
        np.testing.assert_array_equal(output, [15,14,15,15,13,15,15,14,15])

    def test_grid_to_grid_laf_3d(self):
        ivalues = np.ones([5,2,2])*15
        elev_gradient = np.zeros([5,2,2])
        laf_gradient = np.ones([5,2,2])*10
        
        #laf gradient varying in time
        laf_gradient2 = np.ones([5,2,2])*10
        for i in range(len(laf_gradient2[:,0,0])):
            laf_gradient2[i,:,:] = 10 * (1 - i/5.0)

        #Testing constant laf_gradient through time
        output = gridpp.full_gradient(self.igrid, self.ogrid, ivalues, elev_gradient, laf_gradient) 
        np.testing.assert_array_equal(output, [[[10,20,10], [10,20,10], [10,20,10]]]*5)

        #Testing varying laf_gradient through time
        output = gridpp.full_gradient(self.igrid, self.ogrid, ivalues, elev_gradient, laf_gradient2) 
        np.testing.assert_array_equal(output, [[[10,20,10], [10,20,10], [10,20,10]], 
                                            [[11,19,11], [11,19,11], [11,19,11]], 
                                            [[12,18,12], [12,18,12], [12,18,12]],
                                            [[13,17,13], [13,17,13], [13,17,13]],
                                            [[14,16,14], [14,16,14], [14,16,14]]])

    def test_grid_to_grid_elev_3d(self):
        ivalues = np.ones([5,2,2])*15
        elev_gradient = - np.ones([5,2,2]) /100.0 
        laf_gradient = np.zeros([5,2,2])

        elev_gradient2 = np.ones([5,2,2])
        for i in range(len(elev_gradient2[:,0,0])):
            elev_gradient2[i,:,:]  = - 0.01 * (i + 1)

        output = gridpp.full_gradient(self.igrid, self.ogrid, ivalues, elev_gradient, laf_gradient)
        np.testing.assert_array_equal(output,  [[[15,14,15], [15,13,15], [15,14,15]]] * 5 )

        output = gridpp.full_gradient(self.igrid, self.ogrid, ivalues, elev_gradient2, laf_gradient)
        np.testing.assert_array_equal(output,  [[[15,14,15], [15,13,15], [15,14,15]], 
                                            [[15,13,15], [15,11,15], [15,13,15]], 
                                            [[15,12,15], [15,9,15], [15,12,15]], 
                                            [[15,11,15], [15,7,15], [15,11,15]], 
                                            [[15,10,15], [15,5,15], [15,10,15]]])

    def test_grid_to_grid_all_3d(self):
        ivalues = np.ones([5,2,2])*15
        elev_gradient = - np.ones([5,2,2]) / 100.0 
        laf_gradient = np.ones([5,2,2]) * 10

        output = gridpp.full_gradient(self.igrid, self.ogrid, ivalues, elev_gradient, laf_gradient)

        np.testing.assert_array_equal(output, [[[10,19,10], [10,18,10], [10,19,10]]] * 5)

        laf_gradient2 = np.ones([5,2,2])*10
        for i in range(len(laf_gradient2[:,0,0])):
            laf_gradient2[i,:,:] = 10 * (1 - i/5.0) 
        elev_gradient2 = np.ones([5,2,2])
        for i in range(len(elev_gradient2[:,0,0])):
            elev_gradient2[i,:,:]  = - 0.01 * (i + 1) 

        output = gridpp.full_gradient(self.igrid, self.ogrid, ivalues, elev_gradient2, laf_gradient2)
        np.testing.assert_array_equal(output,  [[[10,19,10], [10,18,10], [10,19,10]], 
                                            [[11,17,11], [11,15,11], [11,17,11]], 
                                            [[12,15,12], [12,12,12], [12,15,12]], 
                                            [[13,13,13], [13,9,13], [13,13,13]], 
                                            [[14,11,14], [14,6,14], [14,11,14]]])

    def test_grid_to_point_laf_3d(self):
        pass

    def test_generic(self):
        lats, lons = np.meshgrid([0, 1, 2], [0, 1])
        elevs = np.reshape(np.arange(6), [2, 3])
        grid = gridpp.Grid(lats, lons, elevs)
        gridt = gridpp.Grid(lats.transpose(), lons.transpose(), elevs.transpose())
        lat = 0.9
        lon = 0.9
        elev = 10
        expected = 10-4

        values = np.zeros([2, 3])
        elev_gradient = np.ones([2, 3])
        laf_gradient = np.ones([2, 3])

        self.compute_different_forms(grid, values, elev_gradient, laf_gradient, lat, lon, elev, expected)
        self.compute_different_forms(gridt, values.transpose(), elev_gradient.transpose(), laf_gradient.transpose(), lat, lon, elev, expected)


    def compute_different_forms(self, grid, values, elev_gradient, laf_gradient, lat, lon, elev, expected, *args):
        T = 3
        P = 5

        lats = np.full([P, 1], lat)
        lons = np.full([P, 1], lon)
        elevs = np.full([P, 1], elev)
        ogrid = gridpp.Grid(lats, lons, elevs)
        opoints = gridpp.Points(lats.flatten(), lons.flatten(), elevs.flatten())

        pexpected2 = np.full([P], expected)
        pexpected3 = np.full([T, P], expected)
        gexpected2 = np.full([P, 1], expected)
        gexpected3 = np.full([T, P, 1], expected)

        values3 = np.zeros([T, values.shape[0], values.shape[1]])
        elev_gradient3 = np.zeros([T, values.shape[0], values.shape[1]])
        laf_gradient3 = np.zeros([T, values.shape[0], values.shape[1]])
        for t in range(T):
            values3[t, ...] = values
            laf_gradient3[t, ...] = laf_gradient
            elev_gradient3[t, ...] = elev_gradient

        # Grid to grid 2D
        output = gridpp.full_gradient(grid, ogrid, values, elev_gradient, laf_gradient, *args)
        np.testing.assert_array_almost_equal(output, gexpected2)

        # Grid to grid 3D
        output = gridpp.full_gradient(grid, ogrid, values3, elev_gradient3, laf_gradient3, *args)
        np.testing.assert_array_almost_equal(output, gexpected3)

        # Grid to points 2D
        output = gridpp.full_gradient(grid, opoints, values, elev_gradient, laf_gradient, *args)
        np.testing.assert_array_almost_equal(output, pexpected2)

        # Grid to points 3D
        output = gridpp.full_gradient(grid, opoints, values3, elev_gradient3, laf_gradient3, *args)
        np.testing.assert_array_almost_equal(output, pexpected3)


if __name__ == '__main__':
    unittest.main()
