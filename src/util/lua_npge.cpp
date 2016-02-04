/*
 * NPG-explorer, Nucleotide PanGenome explorer
 * Copyright (C) 2012-2016 Boris Nagaev
 *
 * See the LICENSE file for terms of use.
 */

#include <string>

#include "lua_npge.hpp"
#include "Exception.hpp"

namespace npge {

void init_lua_npge(lua_State* L) {
    std::string code = "npge = require 'npge';"; // FIXME
    int has_error = luaL_dostring(L, code.c_str());
    if (has_error) {
        throw Exception("Can't require module 'npge'.");
    }
}

}
