option(BUILD_CLIENT "Build gridpp command-line client" OFF)
if(BUILD_CLIENT)
    find_package(NetCDF REQUIRED)
    find_package(GSL REQUIRED)

    option(ENABLE_TESTS "build test suite" OFF)
    if (NOT GTEST_DIR)
       message(STATUS "GTEST_DIR is not set, assuming ${CMAKE_CURRENT_SOURCE_DIR}/gtest/")
       set(GTEST_DIR "${CMAKE_CURRENT_SOURCE_DIR}/gtest/" CACHE STRING
          "Set the path to gtest tarball/source directory"
          FORCE)
    endif()
    if(NOT IS_ABSOLUTE "${GTEST_DIR}")
       if (CMAKE_VERSION VERSION_LESS 3.4.3)
          message(FATAL_ERROR "GTEST_DIR path should be absolute")
       else()
          get_filename_component(GTEST_DIR "${GTEST_DIR}" ABSOLUTE
             BASE_DIR ${CMAKE_BINARY_DIR})
       endif()
    endif()

    function(list_sources sources)
       file(GLOB CORESRC "*.cpp")
       file(GLOB CALSRC  "Calibrator/*.cpp")
       file(GLOB DOWNSRC "Downscaler/*.cpp")
       file(GLOB FILESRC "File/*.cpp")
       file(GLOB PARSRC  "ParameterFile/*.cpp")

       set(${sources}
          ${CORESRC}
          ${CALSRC}
          ${SOURCES}
          ${FILESRC}
          ${DOWNSRC}
          ${PARSRC} PARENT_SCOPE)
    endfunction(list_sources)

    function(add_gridpp_library target lib_type)
       add_library(${target} ${lib_type} ${SRC})
       target_include_directories(${target} PUBLIC "${Boost_INCLUDE_DIR}")
       target_include_directories(${target} PUBLIC "${NETCDF_INCLUDE_DIRS}")
       target_include_directories(${target} PUBLIC "${GSL_INCLUDE_DIRS}")
       target_include_directories(${target} PUBLIC "${ARMADILLO_INCLUDE_DIRS}")
       set_target_properties(${target} PROPERTIES OUTPUT_NAME gridpp_client)
    endfunction()

    list_sources(SRC)
    set(DRVSRC  "Driver/Gridpp.cpp")

    add_gridpp_library(gridpp-static STATIC)
    set(GRIDPP_CLIENT_EXE gridpp-client)
    add_executable(${GRIDPP_CLIENT_EXE} ${DRVSRC})
    set_target_properties(${GRIDPP_CLIENT_EXE} PROPERTIES OUTPUT_NAME gridpp)
    set_target_properties(${GRIDPP_CLIENT_EXE} PROPERTIES RUNTIME_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}")
    target_link_libraries(${GRIDPP_CLIENT_EXE} gridpp-static)
    target_link_libraries(${GRIDPP_CLIENT_EXE} "${Boost_LIBRARIES}")
    target_link_libraries(${GRIDPP_CLIENT_EXE} "${NETCDF_LIBRARIES}")
    target_link_libraries(${GRIDPP_CLIENT_EXE} "${GSL_LIBRARIES}")
    target_link_libraries(${GRIDPP_CLIENT_EXE} "${ARMADILLO_LIBRARIES}")
    if (ENABLE_TESTS)
        include(Testing)
        if(CMAKE_PROJECT_NAME STREQUAL PROJECT_NAME AND BUILD_TESTING)
            # Copy the test programs to be next to the gridpp.py package
            # add_custom_target(tests-data
            #     COMMAND ${CMAKE_COMMAND} -E copy_directory
            #     ${CMAKE_CURRENT_SOURCE_DIR}/tests
            #     ${CMAKE_BINARY_DIR}/swig/python/tests)

            # add_custom_target(tests COMMAND nosetests3 swig/python -e memory)
            # add_dependencies(tests tests-data build-python)
            add_custom_target(coverage BYPRODUCTS coverage.info
                COMMAND lcov --directory ${CMAKE_BINARY_DIR}/CMakeFiles/gridpp.dir/ --capture --output-file ${CMAKE_BINARY_DIR}/coverage.info
                COMMAND lcov --remove coverage.info 'client/Testing/*' '/usr/*' --output-file coverage.info
                COMMAND lcov --list coverage.info
                # COMMAND genhtml coverage.info -o coverage
            )
            add_dependencies(coverage tests)
            add_custom_target(coverage-html
                COMMAND genhtml coverage.info -o coverage
            )
            add_dependencies(coverage-html coverage)
        endif()
    endif()
    install(TARGETS gridpp-client RUNTIME DESTINATION bin)
endif()
