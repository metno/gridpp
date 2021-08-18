find_package(PythonLibs REQUIRED)
include_directories(${PYTHON_INCLUDE_PATH})
# swig_add_library(qwegridpp TYPE SHARED LANGUAGE python SOURCES ${CMAKE_SOURCE_DIR}/SWIG/${PROJECT_NAME}.i)

# Files to install with Python
file(GLOB CPP_INCLUDE_FILES ${CMAKE_SOURCE_DIR}/include/*.h)
file(COPY ${CMAKE_CURRENT_SOURCE_DIR}/MANIFEST.in DESTINATION ${CMAKE_CURRENT_BINARY_DIR})
# Configure setup.py and copy to output directory
set(SETUP_PY_IN ${CMAKE_CURRENT_SOURCE_DIR}/setup.in.py)
set(SETUP_PY_OUT ${CMAKE_CURRENT_BINARY_DIR}/setup.py)
configure_file(${SETUP_PY_IN} ${SETUP_PY_OUT})
     #file(COPY ${CMAKE_CURRENT_SOURCE_DIR}/setup.py
     #        DESTINATION ${CMAKE_CURRENT_BINARY_DIR})

# Install target to call setup.py
file(MAKE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/src/api)
file(MAKE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/include)

set_source_files_properties(../SWIG/${PROJECT_NAME}.i PROPERTIES CPLUSPLUS ON)
swig_add_library(gridpp2 TYPE SHARED LANGUAGE python SOURCES ../SWIG/${PROJECT_NAME}.i)
set_target_properties(_gridpp2 PROPERTIES LINKER_LANGUAGE CXX)
add_custom_target(package-python
    COMMAND ${CMAKE_COMMAND} -E copy ${CMAKE_SOURCE_DIR}/README.md ${CMAKE_CURRENT_BINARY_DIR}
    COMMAND ${CMAKE_COMMAND} -E copy ${SOURCES} ${CMAKE_CURRENT_BINARY_DIR}/src
    COMMAND ${CMAKE_COMMAND} -E copy ${CPP_INCLUDE_FILES} ${CMAKE_CURRENT_BINARY_DIR}/include)
add_dependencies(package-python _gridpp2)
