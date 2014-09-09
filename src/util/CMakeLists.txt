file(GLOB sources "*.cpp")

include(CheckCXXSourceCompiles)

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

