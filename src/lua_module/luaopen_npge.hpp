/*
 * NPG-explorer, Nucleotide PanGenome explorer
 * Copyright (C) 2014 Boris Nagaev
 *
 * See the LICENSE file for terms of use.
 */

#ifndef NPGE_LUAOPEN_HPP_
#define NPGE_LUAOPEN_HPP_

#include <lua.hpp>

/** Initialize whole NPGe members in Lua */
extern "C" int luaopen_npge(lua_State* L);

#endif

