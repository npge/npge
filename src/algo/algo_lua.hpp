/*
 * NPG-explorer, Nucleotide PanGenome explorer
 * Copyright (C) 2014 Boris Nagaev
 *
 * See the LICENSE file for terms of use.
 */

#ifndef NPGE_ALGO_LUA_HPP_
#define NPGE_ALGO_LUA_HPP_

#include "lua.hpp"
#include <luabind/luabind.hpp>

#include "global.hpp"

namespace luabind {

typedef std::vector<npge::Processor*> Processors;

template <>
struct default_converter<Processors> :
        native_converter_base<Processors> {
    static int compute_score(lua_State* L, int index);
    Processors from(lua_State* L, int index);
    void to(lua_State* L, const Processors& strings);
};

template <>
struct default_converter<const Processors&> :
        default_converter<Processors> {
};

}

/** Initialize algo/ members in Lua */
extern "C" int init_algo_lua(lua_State* L);

#endif

