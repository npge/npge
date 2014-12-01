file(GLOB sources "*.cpp")
list(APPEND sources
    ${CMAKE_CURRENT_BINARY_DIR}/npge_version.cpp)

include(CheckCXXSourceCompiles)

set(CMAKE_REQUIRED_INCLUDES
    ${LUAIMPL_INCLUDE} ${LUABIND_INCLUDE_DIRS})

CHECK_CXX_SOURCE_COMPILES("#define LUA_COMPAT_ALL\n
    #include <luabind/luabind.hpp>\n
    namespace luabind {
    template <> struct default_converter<int64_t> { }; }
    int main() { }"
    LUABIND_INT64_MISSING)

if (LUABIND_INT64_MISSING)
    add_definitions(-DLUABIND_INT64_MISSING)
endif ()

add_library(npge_util STATIC ${sources})

include_directories(${CMAKE_CURRENT_SOURCE_DIR})
include_directories(${CMAKE_CURRENT_BINARY_DIR})

configure_file(${CMAKE_CURRENT_SOURCE_DIR}/npge_version.cpp.in
    ${CMAKE_CURRENT_BINARY_DIR}/npge_version.cpp)

