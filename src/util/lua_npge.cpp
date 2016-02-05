/*
 * NPG-explorer, Nucleotide PanGenome explorer
 * Copyright (C) 2012-2016 Boris Nagaev
 *
 * See the LICENSE file for terms of use.
 */

#include <iostream>

#include <lua.hpp>

#include "lua_npge.hpp"
#include "throw_assert.hpp"

extern "C" {
int luaopen_npge_cpp(lua_State *L);
}

namespace npge {

void init_lua_npge(lua_State* L) {
    // add npge.cpp module
    lua_getglobal(L, "package");
    ASSERT_EQ(lua_type(L, -1), LUA_TTABLE);
    lua_getfield(L, -1, "loaded");
    luaopen_npge_cpp(L);
    lua_setfield(L, -2, "npge.cpp"); // package.loaded["npge.cpp"] = ...
    lua_pop(L, 2); // package and package.loaded
    // add lua modules of lua-npge
#include "lua_npge.h"
}

}
