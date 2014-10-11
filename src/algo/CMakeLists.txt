include_directories(../util/ ../model/)

configure_file(${CMAKE_CURRENT_SOURCE_DIR}/opts_lib.cpp.in
    ${CMAKE_CURRENT_BINARY_DIR}/opts_lib.cpp)

file(GLOB sources "*.cpp")
list(APPEND sources ${CMAKE_CURRENT_BINARY_DIR}/opts_lib.cpp)

include_directories(${CMAKE_CURRENT_SOURCE_DIR})
include_directories(${CMAKE_CURRENT_BINARY_DIR})

if (USE_LUAJIT)
    add_custom_command(
        OUTPUT ${CMAKE_CURRENT_BINARY_DIR}/lua_lib.h
        COMMAND luajit ${PROJECT_SOURCE_DIR}/src/bin2c.lua
        ${CMAKE_CURRENT_SOURCE_DIR}/lua_lib.lua
            > ${CMAKE_CURRENT_BINARY_DIR}/lua_lib.h
            MAIN_DEPENDENCY
                ${CMAKE_CURRENT_SOURCE_DIR}/lua_lib.lua)
else ()
    add_custom_command(
        OUTPUT ${CMAKE_CURRENT_BINARY_DIR}/lua_lib.h
        COMMAND lua ${PROJECT_SOURCE_DIR}/src/bin2c.lua
        ${CMAKE_CURRENT_SOURCE_DIR}/lua_lib.lua
            > ${CMAKE_CURRENT_BINARY_DIR}/lua_lib.h
            MAIN_DEPENDENCY
                ${CMAKE_CURRENT_SOURCE_DIR}/lua_lib.lua)
endif ()

add_library(npge_algo STATIC ${sources})

set_property(SOURCE lua_lib.cpp APPEND PROPERTY OBJECT_DEPENDS
    ${CMAKE_CURRENT_BINARY_DIR}/lua_lib.h)

