# Author: Shao Zhengjiang<shaozhengjiang@gmail.com>
# Date  : 2017-03-23

find_package(SWIG REQUIRED)
include(${SWIG_USE_FILE})

include_directories(${PROJECT_SOURCE_DIR}/include)

# Read config file.
if (CONFIG)
    file(READ "${PROJECT_SOURCE_DIR}/${CONFIG}" config)
    string(REGEX MATCH "^PYINCLUDE = (.*)\nPYLIB = (.*)\n+$" path ${config})

    # Read PYINCLUDE.
    if (CMAKE_MATCH_1)
        set(PYINCLUDE ${CMAKE_MATCH_1})
    else()
        message(FATAL_ERROR "Missing PYINCLUDE in config file")
    endif()

    # Read PYLIB.
    if (CMAKE_MATCH_2)
        set(PYLIB ${CMAKE_MATCH_2})
    else()
        message(FATAL_ERROR "Missing PYLIB in config file")
    endif()
endif()

# Load PYTHON_INCLUDE_PATH & PYTHON_LIBRARY.
if (PYINCLUDE)
    # Python header file include.
    message( STATUS "Customized Python include path detected" )
    set(PYTHON_INCLUDE_PATH ${PYINCLUDE})
    include_directories(${PYTHON_INCLUDE_PATH})

    # Python library.
    if (NOT PYLIB)
        message(FATAL_ERROR "No Python library provided, use `cmake -DPYLIB=<python-library>`")
    endif()

    message( STATUS "Customized Python library detected" )
    set(PYTHON_LIBRARY ${PYLIB})
else()
    find_package(PythonLibs REQUIRED)
endif()

message(STATUS "Using Python header files in ${PYTHON_INCLUDE_PATH}")
message(STATUS "Using Python library ${PYTHON_LIBRARY}")
# Set flags for swig.
set(CMAKE_SWIG_FLAGS "")

# Set the properties for the interface file.
set_source_files_properties(kmc_functions.i PROPERTIES CPLUSPLUS ON)
set_source_files_properties(kmc_functions.i PROPERTIES SWIG_FLAGS "")

# Add specific compiler flags for macOS
#if (${CMAKE_SYSTEM_NAME} MATCHES "Darwin")
#    add_definitions("-dynamiclib")
#endif()

# Add swig target.
swig_add_module(kmc_functions python kmc_functions.i)

# Link
message(STATUS "Creating makefiles for system: ${CMAKE_SYSTEM}")
swig_link_libraries( kmc_functions ${PYTHON_LIBRARY} src)

