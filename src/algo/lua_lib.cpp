/*
 * NPG-explorer, Nucleotide PanGenome explorer
 * Copyright (C) 2014 Boris Nagaev
 *
 * See the LICENSE file for terms of use.
 */

#include <iostream>

#include <lua.hpp>

#include "lua_lib.hpp"
#include "Meta.hpp"

namespace npge {

void add_lua_lib(Meta* meta) {
    lua_State* L = meta->L();
#include "lua_lib.h"
}

}

