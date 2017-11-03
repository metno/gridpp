Gridded post-processor
======================

.. image:: https://travis-ci.org/metno/gridpp.svg?branch=master
    :target: https://travis-ci.org/metno/gridpp

.. image:: https://coveralls.io/repos/metno/gridpp/badge.svg?branch=master&service=github
    :target: https://coveralls.io/github/metno/gridpp?branch=master 

The program post-processes NetCDF files used at MET-Norway by using various
downscaling and calibration methods. Post-processed forecasts are placed in a
second Netcdf file, which has the desired output grid.

To convert between the two grids, a downscaling method is used. Currently
implemented methods are:

* nearest neighbour

* elevation gradient (interpolation to new elevations using gradients)

* smart neighbours (nearest grid points at the same elevation)

* pressure (interpolation to new elevations using a standard atmosphere)

Calibrators include:

* spatial smoothing (average within a neighbourhood)

* quantile-quantile mapping

* linear regression

* Kriging of biases at points onto a grid (additive and multiplicative)

* ensemble calibration using zero-adjusted Gamma distribution (e.g. for precipitation)

* ensemble calibration using Box-Cox t-distribution (e.g. for windspeed)

* ensemble calibration using Gaussian distribution (e.g. for temperature)

* calculation of precipitation phase, using wetbulb temperature

* calculation of QNH from surface pressure



Installation Instructions
-------------------------

**From source**

1. Ensure the following libraries are installed:

   * Boost libraries
   * Netcdf c++ library
   * libgsl0
   * libblas
   * cmake (version 3 or above)
   * (Optional) Google test library (if developing new code)

2. Download the source code from a release: https://github.com/metno/gridpp/releases

3. Create build directory at the root of the repository

.. code-block:: bash

   mkdir build

4. Run cmake to set up installation

.. code-block:: bash

   cd build
   cmake ..

This will set up gridpp to be installed in /usr/local/bin on Ubuntu. To specify a custom installation path, use:

.. code-block:: bash

   cmake .. -DCMAKE_INSTALL_PREFIX=<custom path>

5. Run cmake to install

.. code-block:: bash

   cmake --build .
   cmake --build . --target install


**From debian packages**

Follow instructions here: https://launchpad.net/~metno/+archive/ubuntu/gridpp


Running the program
-------------------
To see program options, run:

.. code-block:: bash

   gridpp

To test the program on a fake dataset, run:

.. code-block:: bash

   gridpp testing/files/10x10.nc testing/files/10x10_copy.nc\
      -v T -d gradient\
      -v Precip -d smart numSmart=3 searchRadius=3
   ncview testing/files/10x10_copy.nc



Running on multiple cores
-------------------------
To run using 8 threads:

.. code-block:: bash

   export OMP_NUM_THREADS=8
   gridpp ...



Minimizing memory usage
-----------------------
Run the program in sequence for each variable:

.. code-block:: bash

   gridpp input output -v T ...
   gridpp input output -v Precip ...
   gridpp input output -v RH ...



Copyright and license
---------------------
Copyright (C) 2015-2017 MET Norway. Gridded post-processor is licensed under `GPL
version 2 <https://github.com/metno/gridpp/blob/master/LICENSE>`_ or (at
your option) any later version.

Contact
-------
| E-mail: `thomasn@met.no <mailto:thomasn@met.no>`_
