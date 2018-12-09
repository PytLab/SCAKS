# Author: Shao Zhengjiang<shaozhengjiang@gmail.com>
# Date  : 2017-03-23

cmake_minimum_required(VERSION 2.8)
project(Kynetix)

set(CMAKE_VERBOSE_MAKEFILE on)

# Set compiler related parameters.
set( CMAKE_CXX_FLAGS "-O3 -fPIC -std=c++0x -Wall -Werror" )

# Add the source target to be build.
add_subdirectory(src)

# Add swig wrap targets to be build.
add_subdirectory(wrap)

# Install the target
install(FILES ${PROJECT_BINARY_DIR}/wrap/kmc_functions.py ${PROJECT_BINARY_DIR}/wrap/_kmc_functions.so
    DESTINATION ${PROJECT_SOURCE_DIR}/../scaks/solvers/kmc_plugins/plugin_backends/)

