/*
 * NPG-explorer, Nucleotide PanGenome explorer
 * Copyright (C) 2014 Boris Nagaev
 *
 * See the LICENSE file for terms of use.
 */

#include <iostream>

#include <lua.hpp>

#include "lua_lib.hpp"
#include "lua_lib.lua"
#include "Meta.hpp"

namespace npge {

void add_lua_lib(Meta* meta) {
    if (luaL_dostring(meta->L(), meta_lua)) {
        // error
        std::cerr << lua_tostring(meta->L(), -1) << "\n";
    }
}

}

