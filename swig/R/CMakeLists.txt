find_package(R REQUIRED)
include_directories(${R_INCLUDE_DIR})

set(CMAKE_SWIG_FLAGS "")
set_source_files_properties(../${PROJECT_NAME}.i PROPERTIES CPLUSPLUS ON)

include_directories(../../include)

# Add swig module
swig_add_library(rgridpp TYPE SHARED LANGUAGE r SOURCES ../${PROJECT_NAME}.i ${SOURCES})
swig_link_libraries(rgridpp ${PROJECT_NAME} ${R_LIBRARIES})
#set_property(SOURCE ../gridpp.i PROPERTY SWIG_MODULE_NAME gridpp)
#set_property(TARGET rgridpp PROPERTY SWIG_USE_TARGET_INCLUDE_DIRECTORIES ON)
set_property(TARGET rgridpp PROPERTY OUTPUT_NAME gridpp)

# Files to install with Python
#set(PYTHON_INSTALL_FILES
#        ${CMAKE_CURRENT_BINARY_DIR}/gridpp.py
#        ${CMAKE_CURRENT_BINARY_DIR}/_gridpp.so)


# Declare install target for python
#install(TARGETS swig_example
#        COMMAND "${PYTHON_EXECUTABLE} setup.py"
#        COMPONENT swig-python)

# Install target to call setup.py
#target_link_libraries(rgridpp
#  PRIVATE
#  gridpp
#  ${R_LIBRARIES}
#  )
add_custom_target(build-r
        DEPENDS rgridpp
        )
