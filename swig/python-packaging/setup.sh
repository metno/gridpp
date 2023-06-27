mkdir build
cd build
cmake .. -DBUILD_PACKAGE=ON
VERBOSE=1 make package-python
