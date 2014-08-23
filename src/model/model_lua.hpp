/*
 * NPG-explorer, Nucleotide PanGenome explorer
 * Copyright (C) 2014 Boris Nagaev
 *
 * See the LICENSE file for terms of use.
 */

#ifndef NPGE_MODEL_LUA_HPP_
#define NPGE_MODEL_LUA_HPP_

#include "lua.hpp"

/** Initialize model/ members in Lua */
extern "C" int init_model_lua(lua_State* L);

#endif

