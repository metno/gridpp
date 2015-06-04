Notes to developers
===================

Adding new methods
------------------
Here is a recipe for adding a new calibration method called MyMethod:

1) Create src/Calibrator/MyMethod.h and src/Calibrator/MyMethod.cpp. You can base the files on one of the other
   calibrators. Make sure the class name is CalibratorMyMethod.
   The header file should have the include guards:

  .. code-block:: c++

    #ifndef CALIBRATOR_MY_METHOD_H
    #define CALIBRATOR_MY_METHOD_H

2) Follow the instructions in src/Calibrator/Calibrator.h on what methods you must implement.
3) Add #include "MyScheme.h" at the bottom of Calibrator.h.
4) Define how your method is instantiated inside the getScheme function in src/Calibrator/Calibrator.cpp:

  .. code-block:: c++

    else if(iName == "myMethod") {
      CalibratorMyMethod* c = new CalibratorMyMethod(Variable::getType(variable));
    }

5) Add the following to src/Driver/Gridpp.cpp so that the description of your method shows up in gridpp's
   usage information:

  .. code-block:: c++

    std::cout << CalibratorMyMethod::description();

Debugging
---------
gridpp can be compiled with debugging symbols by executing make gridpp_debug. The executable gridpp_debug can
then be used by a debugger such as gdb.

Parallelization
---------------
Add openMP directives where it makes sense to parallelize. Ensure that any data retrival calls (e.g.
calls to File::getField) are made before the directives otherwise a runtime error occurs.

Testing of code
---------------
All test code is placed in src/Testing. Use one file for each class. To test your newly added class
create a test file called src/Testing/CalibratorMyMethod.cpp. Follow the format in the other test files.
Aim to test your class thoroughly. There is a test file located in testing/files/10x10.nc. You can test
that your calibration produces results as expected on data in this file.

To check how much of the code is tested you can run the following in
the top of the repository:

  .. code-block:: bash

    make coverage

Open coverage/index.html in the browser to see the coverage results.

Releasing a new version
-----------------------
1) Commit your changes.
2) Determine the new tag. Use "git tag" to find the latest tag and increment according to
   http://semver.org/. The tag should be in the form v0.1.1.
3) Edit src/Version.h to indicate the new version (without the v-prefix). This value is used
   by gridpp --version.
4) Update the debian changelog by putting in the version number and filling in a description
       of the new release.

  .. code-block:: bash

     dch -i

5) Run the test suite and make sure there are no failures. The tests should be rerun here
   in case editing srs/Version.h caused compile errors.

  .. code-block:: bash

    make test

6) Commit the release information

  .. code-block:: bash

    git commit debian/changelog src/Version.h

7) Tag the version in git (using the previously determined tag)

  .. code-block:: bash

     git tag <tag including the v-prefix>

7) Push the release to the repository

  .. code-block:: bash

     git push --tags origin master
