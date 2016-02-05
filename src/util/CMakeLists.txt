file(GLOB sources "*.cpp")
list(APPEND sources
    ${CMAKE_CURRENT_BINARY_DIR}/npge_version.cpp)

file(GLOB lua_npge_sources "../lua-npge/src/npge/cpp/*.cpp")

set(NPGE_LUA_CMD "lua" CACHE STRING "Lua command")

execute_process(
    COMMAND ${NPGE_LUA_CMD} -e x=x
    RESULT_VARIABLE HAS_LUA_CMD
)
if (NOT (HAS_LUA_CMD EQUAL 0))
    message(SEND_ERROR
            "Specify -DNPGE_LUA_CMD=/path/to/lua.exe")
endif()

add_custom_command(
    OUTPUT ${CMAKE_CURRENT_BINARY_DIR}/lua_npge.lua
    COMMAND ${NPGE_LUA_CMD} ${PROJECT_SOURCE_DIR}/src/statify_lua-npge.lua
    ${PROJECT_SOURCE_DIR}
        > ${CMAKE_CURRENT_BINARY_DIR}/lua_npge.lua
)

add_custom_command(
    OUTPUT ${CMAKE_CURRENT_BINARY_DIR}/lua_npge.h
    COMMAND ${NPGE_LUA_CMD} ${PROJECT_SOURCE_DIR}/src/bin2c.lua
    ${CMAKE_CURRENT_BINARY_DIR}/lua_npge.lua
        > ${CMAKE_CURRENT_BINARY_DIR}/lua_npge.h
    MAIN_DEPENDENCY
    ${CMAKE_CURRENT_BINARY_DIR}/lua_npge.lua
)

set_property(SOURCE lua_npge.cpp APPEND PROPERTY OBJECT_DEPENDS
    ${CMAKE_CURRENT_BINARY_DIR}/lua_npge.h
)

include(CheckCXXSourceCompiles)

set(CMAKE_REQUIRED_INCLUDES
    ${LUAIMPL_INCLUDE} ${LUABIND_INCLUDE_DIRS})

CHECK_CXX_SOURCE_COMPILES("#define LUA_COMPAT_ALL\n
    #include <stdint.h>\n
    #include <luabind/luabind.hpp>\n
    namespace luabind {
    template <> struct default_converter<int64_t> { }; }
    int main() { }"
    LUABIND_INT64_MISSING)

if (LUABIND_INT64_MISSING)
    add_definitions(-DLUABIND_INT64_MISSING)
endif ()

add_library(npge_util STATIC ${sources} ${lua_npge_sources})

include_directories(${CMAKE_CURRENT_SOURCE_DIR})
include_directories(${CMAKE_CURRENT_BINARY_DIR})

configure_file(${CMAKE_CURRENT_SOURCE_DIR}/npge_version.cpp.in
    ${CMAKE_CURRENT_BINARY_DIR}/npge_version.cpp)

