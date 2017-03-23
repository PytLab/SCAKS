include_directories(${PROJECT_SOURCE_DIR}/include)

file(GLOB Cppfiles *.cpp)
add_library(src  ${Cppfiles})

