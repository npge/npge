include_directories(../util/ ../model/)

configure_file(${CMAKE_CURRENT_SOURCE_DIR}/opts_lib.cpp.in
    ${CMAKE_CURRENT_BINARY_DIR}/opts_lib.cpp)

file(GLOB sources "*.cpp")
list(APPEND sources ${CMAKE_CURRENT_BINARY_DIR}/opts_lib.cpp)

include_directories(${CMAKE_CURRENT_SOURCE_DIR})

add_library(npge_algo STATIC ${sources})

