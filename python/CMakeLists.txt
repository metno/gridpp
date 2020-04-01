find_package(PythonLibs REQUIRED)
include_directories(${PYTHON_INCLUDE_PATH})
# swig_add_library(qwegridpp TYPE SHARED LANGUAGE python SOURCES ${CMAKE_SOURCE_DIR}/SWIG/${PROJECT_NAME}.i)

# Files to install with Python
file(GLOB CPP_SOURCE_FILES ${CMAKE_SOURCE_DIR}/src/*.cpp ${CMAKE_SOURCE_DIR}/src/*.c)
file(GLOB CPP_INCLUDE_FILES ${CMAKE_SOURCE_DIR}/include/*.h)
file(COPY ${CMAKE_CURRENT_SOURCE_DIR}/MANIFEST.in DESTINATION ${CMAKE_CURRENT_BINARY_DIR})
#file(COPY ${CPP_SOURCE_FILES} DESTINATION ${CMAKE_CURRENT_BINARY_DIR}/src)
#file(COPY ${CPP_INCLUDE_FILES} DESTINATION ${CMAKE_CURRENT_BINARY_DIR}/include)
#file(COPY ${SWIG_INTERFACE} DESTINATION ${CMAKE_CURRENT_BINARY_DIR})
# Configure setup.py and copy to output directory
set(SETUP_PY_IN ${CMAKE_CURRENT_SOURCE_DIR}/setup.in.py)
set(SETUP_PY_OUT ${CMAKE_CURRENT_BINARY_DIR}/setup.py)
configure_file(${SETUP_PY_IN} ${SETUP_PY_OUT})
     #file(COPY ${CMAKE_CURRENT_SOURCE_DIR}/setup.py
     #        DESTINATION ${CMAKE_CURRENT_BINARY_DIR})

# Install target to call setup.py
MESSAGE(${CPP_SOURCE_FILES})
file(MAKE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/src/api)
file(MAKE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/include)

#configure_file(${SWIG_INTERFACE} ${CMAKE_CURRENT_BINARY_DIR}/)
#foreach(FILE ${CPP_SOURCE_FILES})
#    get_filename_component(F ${FILE} NAME)
#    # MESSAGE(${F})
#    # MESSAGE(${CMAKE_CURRENT_BINARY_DIR})
#    # MESSAGE(${CMAKE_SOURCE_DIR}/src/${F})
#    # MESSAGE(${CMAKE_CURRENT_BINARY_DIR}/src/${F})
#    MESSAGE("COPYING ${F}")
#    configure_file(${CMAKE_SOURCE_DIR}/src/${F} ${CMAKE_CURRENT_BINARY_DIR}/src/${F})
#endforeach()
#foreach(FILE ${CPP_INCLUDE_FILES})
#    get_filename_component(F ${FILE} NAME)
#    configure_file(${CMAKE_SOURCE_DIR}/include/${F} ${CMAKE_CURRENT_BINARY_DIR}/include/${F})
#endforeach()

add_custom_target(package-python
    COMMAND ${CMAKE_COMMAND} -E copy ${CMAKE_SOURCE_DIR}/README.rst ${CMAKE_CURRENT_BINARY_DIR}
    COMMAND ${CMAKE_COMMAND} -E copy ${CMAKE_BINARY_DIR}/SWIG/python/gridppPYTHON_wrap.cxx ${CMAKE_CURRENT_BINARY_DIR}/src
    COMMAND ${CMAKE_COMMAND} -E copy ${CMAKE_BINARY_DIR}/SWIG/python/gridpp.py ${CMAKE_CURRENT_BINARY_DIR}
    COMMAND ${CMAKE_COMMAND} -E copy ${APISRC} ${CMAKE_CURRENT_BINARY_DIR}/src/api
    COMMAND ${CMAKE_COMMAND} -E copy ${CPP_INCLUDE_FILES} ${CMAKE_CURRENT_BINARY_DIR}/include
    COMMAND ${CMAKE_COMMAND} -E copy ${SWIG_INTERFACE} ${CMAKE_CURRENT_BINARY_DIR}
    COMMAND python3 ${SETUP_PY_OUT} sdist)
add_dependencies(package-python build-python)
