include_directories(../util/)

file(GLOB sources "*.cpp")

add_library(npge_model STATIC ${sources})

