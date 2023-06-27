cmake_minimum_required(VERSION 3.0)

# Determine gridpp version
file (READ "${CMAKE_CURRENT_SOURCE_DIR}/include/gridpp.h" VERSION0)
project(gridpp)
string(REGEX REPLACE "^.+define GRIDPP_VERSION \"([^\n]+)\".+" "\\1" PROJECT_VERSION ${VERSION0})

set(CMAKE_MODULE_PATH ${CMAKE_CURRENT_SOURCE_DIR}/cmake)

option(BUILD_PACKAGE "Build package" OFF)

if (BUILD_PACKAGE)

else()
    # Boost 1.75 and above requires c++14 for the geometry library. For creating the python packages on
    # CentOS5, c++11 is more readily available, so this is preferred.
    find_package(Boost REQUIRED)
    find_package(Armadillo 6.5 REQUIRED)
    if(${Boost_MINOR_VERSION} GREATER 74)
        set(CMAKE_CXX_STANDARD 14)
    else()
        set(CMAKE_CXX_STANDARD 11)
    endif()

    # Set compiler flags
    enable_language(CXX)

    set(CMAKE_CXX_STANDARD_REQUIRED ON)
    set(CMAKE_CXX_EXTENSIONS OFF)
    if (NOT CMAKE_BUILD_TYPE)
       message(STATUS "CMAKE_BUILD_TYPE is not set, switching to RELEASE")
       set(CMAKE_BUILD_TYPE "RELEASE" CACHE STRING
          "Choose the type of build, options are: NONE DEBUG RELEASE."
          FORCE)
    endif()
    string(TOUPPER ${CMAKE_BUILD_TYPE} CMAKE_BUILD_TYPE)

    if(APPLE)
        set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -Wfatal-errors -Wno-reorder -Wno-unused-variable -Wno-unused-value -Wno-narrowing -Wno-sign-compare -Wno-Wdelete-non-virtual-dtor -Wno-c++11-narrowing")
    else()
        # SWIG generates a lot of warnings
        set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -Wfatal-errors -Wno-reorder -Wno-unused-variable -Wno-unused-value -Wno-narrowing -Wno-sign-compare -Wno-delete-non-virtual-dtor -Wl,-init,initialize_omp")
    endif()

    if (CMAKE_BUILD_TYPE MATCHES "DEBUG")
       set(CMAKE_CXX_FLAGS_DEBUG "-O0 -fPIC -g -pg -rdynamic --coverage")
       set(CMAKE_LFLAGS_DEBUG " --coverage")
       add_definitions("-DDEBUG")
    elseif(CMAKE_BUILD_TYPE MATCHES "RELEASE")
       find_package(OpenMP REQUIRED)
       if(APPLE)
            set(CMAKE_CXX_FLAGS_RELEASE "-O3 -fPIC")
        else()
            set(CMAKE_CXX_FLAGS_RELEASE "-O3 -fPIC ${OpenMP_CXX_FLAGS}")
        endif()
    else()
        message(FATAL_ERROR "Unknown release type")
    endif()

    # Required packages

    add_library(gridpp ${SOURCES})
    target_link_libraries(gridpp "${Boost_LIBRARIES}")
    target_link_libraries(gridpp "${ARMADILLO_LIBRARIES}")
    target_include_directories(gridpp PUBLIC "${Boost_INCLUDE_DIR}")
    target_include_directories(gridpp PUBLIC "${ARMADILLO_INCLUDE_DIRS}")

endif()
include_directories(./include)
file(GLOB SOURCES "src/api/*.cpp")

# Documentation
option(BUILD_DOC "Build documentation" ON)
if (BUILD_DOC)
    find_package(Doxygen)
    if (DOXYGEN_FOUND)
        set(DOXYGEN_IN ${CMAKE_CURRENT_SOURCE_DIR}/docs/config.in)
        set(DOXYGEN_OUT ${CMAKE_CURRENT_BINARY_DIR}/Doxyfile)
        configure_file(${DOXYGEN_IN} ${DOXYGEN_OUT} @ONLY)
        configure_file("${CMAKE_CURRENT_SOURCE_DIR}/docs/DoxygenLayout.xml" "DoxygenLayout.xml" @ONLY)
        add_custom_target(docs ALL
            COMMAND ${DOXYGEN_EXECUTABLE} ${DOXYGEN_OUT}
            WORKING_DIRECTORY
            ${CMAKE_CURRENT_BINARY_DIR}
            COMMENT "Generating API
            documentation with Doxygen"
            VERBATIM )
    else()
        message(STATUS "Documentation generation disabled. Doxygen needs to be installed to generate documentation.")
    endif()
endif()

add_subdirectory(src/client)
file(GLOB SWIG_INTERFACE "swig/gridpp.i")
add_subdirectory(swig)
