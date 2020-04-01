find_package(PythonLibs REQUIRED)
include_directories(${PYTHON_INCLUDE_PATH})

set(CMAKE_SWIG_FLAGS "")
set_source_files_properties(../${PROJECT_NAME}.i PROPERTIES CPLUSPLUS ON)
include_directories(../../include)

# Add swig module
#set(CMAKE_SWIG_FLAGS -interface _pygridpp)
swig_add_library(gridpp TYPE SHARED LANGUAGE python SOURCES ../${PROJECT_NAME}.i)
swig_link_libraries(gridpp ${PROJECT_NAME} ${PYTHON_LIBRARIES})
#set_property(TARGET _pygridpp PROPERTY OUTPUT_NAME gridpp)

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
        COMMAND python3 ${SETUP_PY_OUT} install)

add_custom_target(install-python-user
        DEPENDS _gridpp
        COMMAND python3 ${SETUP_PY_OUT} install --user)

add_custom_target(build-python
        DEPENDS _gridpp
        COMMAND python3 ${SETUP_PY_OUT} build)

add_custom_target(develop-python
        DEPENDS _gridpp
        COMMAND python3 ${SETUP_PY_OUT} develop)

add_custom_target(develop-python-user
        DEPENDS _gridpp
        COMMAND python3 ${SETUP_PY_OUT} develop --user)

