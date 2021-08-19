find_package(PythonLibs REQUIRED)
find_package(PythonInterp REQUIRED)
include_directories(${PYTHON_INCLUDE_PATH})

if(SWIG_VERSION GREATER_EQUAL 4)
    set(CMAKE_SWIG_FLAGS "-doxygen")
endif()
set_source_files_properties(../${PROJECT_NAME}.i PROPERTIES CPLUSPLUS ON)
include_directories(../../include)

# Get Numpy include directory
# http://mad.web.cern.ch/mad/releases/madng/madng-git/lib/nlopt-2.6.1/cmake/FindNumPy.cmake
execute_process(COMMAND "${PYTHON_EXECUTABLE}" "-c"
    "import numpy as np; print(np.__version__); print(np.get_include());"
    RESULT_VARIABLE _NUMPY_SEARCH_SUCCESS
    OUTPUT_VARIABLE _NUMPY_VALUES_OUTPUT
    ERROR_VARIABLE _NUMPY_ERROR_VALUE
    OUTPUT_STRIP_TRAILING_WHITESPACE)

if(_NUMPY_SEARCH_SUCCESS MATCHES 0)
  # Convert the process output into a list
  string(REGEX REPLACE ";" "\\\\;" _NUMPY_VALUES ${_NUMPY_VALUES_OUTPUT})
  string(REGEX REPLACE "\n" ";" _NUMPY_VALUES ${_NUMPY_VALUES})
  list(GET _NUMPY_VALUES 0 NUMPY_VERSION)
  list(GET _NUMPY_VALUES 1 NUMPY_INCLUDE_DIRS)

  string(REGEX MATCH "^[0-9]+\\.[0-9]+\\.[0-9]+" _VER_CHECK "${NUMPY_VERSION}")
  if("${_VER_CHECK}" STREQUAL "")
      # The output from Python was unexpected. Raise an error always
      # here, because we found NumPy, but it appears to be corrupted somehow.
      message(FATAL_ERROR
          "Requested version and include path from NumPy, got instead:\n${_NUMPY_VALUES_OUTPUT}\n")
      return()
  endif()

  # Make sure all directory separators are '/'
  string(REGEX REPLACE "\\\\" "/" NUMPY_INCLUDE_DIRS ${NUMPY_INCLUDE_DIRS})

  # Get the major and minor version numbers
  string(REGEX REPLACE "\\." ";" _NUMPY_VERSION_LIST ${NUMPY_VERSION})
  list(GET _NUMPY_VERSION_LIST 0 NUMPY_VERSION_MAJOR)
  list(GET _NUMPY_VERSION_LIST 1 NUMPY_VERSION_MINOR)
  list(GET _NUMPY_VERSION_LIST 2 NUMPY_VERSION_PATCH)
  string(REGEX MATCH "[0-9]*" NUMPY_VERSION_PATCH ${NUMPY_VERSION_PATCH})
  math(EXPR NUMPY_VERSION_DECIMAL
      "(${NUMPY_VERSION_MAJOR} * 10000) + (${NUMPY_VERSION_MINOR} * 100) + ${NUMPY_VERSION_PATCH}")
endif()

# handle REQUIRED and QUIET options
include (FindPackageHandleStandardArgs)
find_package_handle_standard_args (NumPy REQUIRED_VARS NUMPY_INCLUDE_DIRS VERSION_VAR NUMPY_VERSION)

# Include Numpy Include directory
include_directories(${NUMPY_INCLUDE_DIRS})

# Add swig module
#set(CMAKE_SWIG_FLAGS -interface _pygridpp)
swig_add_library(gridpp TYPE SHARED LANGUAGE python SOURCES ../${PROJECT_NAME}.i)
#swig_link_libraries(gridpp ${PROJECT_NAME} ${PYTHON_LIBRARIES})
swig_link_libraries(gridpp gridpp ${PYTHON_LIBRARIES})
#set_property(TARGET _pygridpp PROPERTY OUTPUT_NAME gridpp)
if(APPLE)
    set_property(TARGET _gridpp PROPERTY SWIG_USE_TARGET_INCLUDE_DIRECTORIES ON)
    set_target_properties(_gridpp PROPERTIES
       SUFFIX ".so"
       INSTALL_RPATH"@loader_path;@loader_path/../../${PROJECT_NAME}/.libs"
       )
     set_property(TARGET _gridpp APPEND PROPERTY
       LINK_FLAGS "-flat_namespace -undefined suppress"
       )
endif()


# Files to install with Python
set(PYTHON_INSTALL_FILES
        ${CMAKE_CURRENT_BINARY_DIR}/gridpp.py
        ${CMAKE_CURRENT_BINARY_DIR}/_gridpp.so)

# Configure setup.py and copy to output directory
set(SETUP_PY_IN ${CMAKE_CURRENT_SOURCE_DIR}/setup.in.py)
set(SETUP_PY_OUT ${CMAKE_CURRENT_BINARY_DIR}/setup.py)
configure_file(${SETUP_PY_IN} ${SETUP_PY_OUT})

#set(INIT_PY_IN ${CMAKE_CURRENT_SOURCE_DIR}/version.in.py)
#set(INIT_PY_OUT ${CMAKE_CURRENT_BINARY_DIR}/gridpp/__init__.py)
#configure_file(${INIT_PY_IN} ${INIT_PY_OUT})

# Install target to call setup.py
add_custom_target(install-python
        DEPENDS _gridpp
        COMMAND ${PYTHON_EXECUTABLE} ${SETUP_PY_OUT} install)

add_custom_target(install-python-user
        DEPENDS _gridpp
        COMMAND ${PYTHON_EXECUTABLE} ${SETUP_PY_OUT} install --user)

add_custom_target(build-python
        DEPENDS _gridpp
        COMMAND ${PYTHON_EXECUTABLE} ${SETUP_PY_OUT} build)

add_custom_target(develop-python
        DEPENDS _gridpp
        COMMAND ${PYTHON_EXECUTABLE} ${SETUP_PY_OUT} develop)

add_custom_target(develop-python-user
        DEPENDS _gridpp
        COMMAND ${PYTHON_EXECUTABLE} ${SETUP_PY_OUT} develop --user)

if (ENABLE_TESTS)
    if(CMAKE_PROJECT_NAME STREQUAL PROJECT_NAME AND BUILD_TESTING)
        # Copy the test programs to be next to the gridpp.py package
        add_custom_target(tests-data
            COMMAND ${CMAKE_COMMAND} -E copy_directory
            ${CMAKE_SOURCE_DIR}/tests
            ${CMAKE_CURRENT_BINARY_DIR}/tests)

        add_custom_target(tests COMMAND nosetests3 -e memory)
        add_dependencies(tests tests-data build-python)
    endif()
endif()
