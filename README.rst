Gridded post-processor
======================

.. image:: https://travis-ci.org/metno/gridpp.svg?branch=master
    :target: https://travis-ci.org/metno/gridpp

The program post-processes an AROME or ECMWF NetCDF file by using various
downscaling and calibration methods. Post-processed forecasts are placed in a
second Netcdf file, which has the desired output grid.

To convert between the two grids, a downscaling method is used. Currently
implemented methods are:

* nearest neighbour

* elevation gradient

* smart neighbours (nearest grid points at the same elevation)

Calibrators include:

* spatial smoothing

* ensemble calibration of precipitation using zero-adjusted Gamma distribution

* calculation of precipitation phase, using wetbulb temperature




Installation Instructions
-------------------------

1. Ensure the following libraries are installed:

   * Boost libraries
   * Netcdf c++ library
   * (Optional) Google test library (if developing new code)

2. Edit CC, CFLAGS_O, IFLAGS, and LFLAGS in makefile

3. Run 'make'



Running the program
-------------------
To see program options, run:

.. code-block:: bash

   ./gridpp

To test the program on a fake dataset, run:

.. code-block:: bash

   ./gridpp testing/files/10x10.nc testing/files/10x10_copy.nc -v T -d gradient -v Precip -d smart numSmart=3 searchRadius=3
   ncview testing/files/10x10_copy.nc



Example
-------
To test the program on a real operational AROME file, follow these steps:

1. Make a copy of a file to post-process:

.. code-block:: bash

   cp /starc/DNMI_AROME_METCOOP/2015/01/01/AROME_MetCoOp_00_DEF.nc_20150101 localcopy.nc

2. Create a smoothed precipitation field in localcopy.nc:

.. code-block:: bash

   ./gridpp /starc/DNMI_AROME_METCOOP/2015/01/01/AROME_MetCoOp_00_DEF.nc_20150101 localcopy.nc -v Precip -c smooth smoothRadius=10



Running on multiple cores
-------------------------
To run using 8 threads:

.. code-block:: bash

   export OMP_NUM_THREADS=8
   ./gridpp ...



Minimizing memory usage
-----------------------
Run the program in sequence for each variable:

.. code-block:: bash

   ./gridpp input output -v T ...
   ./gridpp input output -v Precip ...
   ./gridpp input output -v RH ...



Adding new code
---------------
* Add a new downscaling scheme by creating a new file in ./src/Downscaler/.
  Inherit from the Downscaler abstract base class (./src/Dowscaler/Downscaler.h)
  and implement any pure virtual functions.
* Instantiate the class in getScheme() in ./src/Downscaler/Downscaler.cpp.
* Use a similar procedure for adding a new calibrator, or new file type.



Testing of code
---------------
* Set DEBUG to 1 in makefile, and run make test
* Execute runAllTests.csh
* Unit tests are placed in src/Testing/, one file for each class that is tested.
* Convenient input files for testing are located in testing/files/
* New downscalers and calibrators should produce valid results for the special
  test file 'testing/files/10x10.nc'.

Copyright and license
---------------------
Copyright (C) 2015 MET Norway. Gridded post-processor is licensed under `GPL
version 2 <https://github.com/metno/gridpp/blob/master/LICENSE>`_ or (at
your option) any later version.

Contact
-------
| MET Norway
| Postboks 43 Blindern
| NO-0313 OSLO
|
| Website: http://met.no/
| E-mail: `post@met.no <mailto:post@met.no>`_
