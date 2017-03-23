# Author: Shao Zhengjiang<shaozhengjiang@gmail.com>
# Date  : 2017-03-23

include_directories(${PROJECT_SOURCE_DIR}/include)

file(GLOB Cppfiles *.cpp)
add_library(src  ${Cppfiles})

