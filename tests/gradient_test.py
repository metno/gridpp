import gridpp
import unittest
import numpy as np 
import random
import matplotlib.pyplot as plt 


def display_gradients(size = 10):
	laf_values = np.zeros([size,size], float)
	t_values = np.zeros([size,size], float)
	for i in range(len(laf_values[:,0])):
		for j in range(len(laf_values[0,:])):
			if i + j < 7:
				laf_values[i,j] = 0
				t_values[i,j] = 10# + random.randrange(-100, 100)/100.0
			elif i + j >= 7:
				laf_values[i,j] = (i + j)/(size**2/5.0)# + random.randrange(-100, 100)/1000.0
				t_values[i,j] = size +  (i + j)/(size/5.0)# + random.randrange(-100, 100)/100.0
			else:
				pass

	gradient = gridpp.calc_gradient(laf_values, t_values, 5)

	fig, (ax1, ax2, ax3) = plt.subplots(figsize=(13, 3), ncols=3)
	laf = ax1.imshow(laf_values, cmap = 'Blues_r', vmin= 0.0, vmax=1.0)
	cbar = fig.colorbar(laf, ax=ax1, extend= 'both')
	cbar.minorticks_on()
	ax1.set_title('Base')

	t = ax2.imshow(t_values, cmap = 'RdBu_r')
	cbar2 = fig.colorbar(t, ax=ax2, extend= 'both')
	cbar2.minorticks_on()
	ax2.set_title('Values')

	grdnt = ax3.imshow(gradient, cmap = 'Greys')
	cbar3 = fig.colorbar(grdnt, ax=ax3, extend= 'both')
	cbar3.minorticks_on()
	ax3.set_title('gradient')

	plt.savefig('output_gradient.png')
	plt.close()
	print(gradient)
	return

class Test(unittest.TestCase):	
	def setup(self):
		""" Setup 10x10 matrix with gradual values with/without randomness"""
		self.size = 10
		self.input_base = np.zeros([self.size,self.size], float)
		self.input_values = np.zeros([self.size,self.size], float)

		""" """
		for i in range(len(self.input_base[:,0])):
			for j in range(len(self.input_base[0,:])):
				if i + j < 7:
					self.input_base[i,j] = 0
					self.input_values[i,j] = 10  #+ random.randrange(-100, 100)/100.0
				elif i + j >= 7:
					self.input_base[i,j] = (i + j)/(self.size**2/5.0)# + random.randrange(-100, 100)/1000.0
					self.input_values[i,j] = self.size + (i + j)/(self.size/5.0)# + random.randrange(-100, 100)/100.0

	def test_size_of_inputs(self):
		""" Checks that the input size of input and output gradient matrix are at the same size""" 
		np.testing.assert_equal(self.input_base.shape, gridpp.calc_gradient(self.input_base, self.input_values, 3).shape)
		np.testing.assert_equal(self.input_values.shape, gridpp.calc_gradient(self.input_base, self.input_values, 3).shape)

		np.testing.assert_array_equal(self.input_base, gridpp.calc_gradient(self.input_base, self.input_values, 3))
		np.testing.assert_array_equal(self.input_values, gridpp.calc_gradient(self.input_base, self.input_values, 3))

	def test_gradient(self):
		""" tests that inputs give the correct gradient"""
		output = gridpp.calc_gradient(self.input_base, self.input_values, 3)
		for i in range(len(self.input_base[0,:])):
			for j in range(len(self.input_base[:,0])):
				np.testing.assert_almost_equal(output[i,j], 10)

	def test_empty_input(self):
		np.testing.assert_array_almost_equal(gridpp.calc_gradient([],[]), [])

if __name__ == '__main__':
	unittest.main()

