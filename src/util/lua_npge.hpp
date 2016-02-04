/*
 * NPG-explorer, Nucleotide PanGenome explorer
 * Copyright (C) 2012-2016 Boris Nagaev
 *
 * See the LICENSE file for terms of use.
 */

#ifndef NPGE_LUA_NPGE_HPP_
#define NPGE_LUA_NPGE_HPP_

#include <lua.hpp>

namespace npge {

/** Add global var "npge" and module "npge" */
void init_lua_npge(lua_State* L);

}

#endif
