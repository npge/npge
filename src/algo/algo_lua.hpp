/*
 * NPG-explorer, Nucleotide PanGenome explorer
 * Copyright (C) 2012-2016 Boris Nagaev
 *
 * See the LICENSE file for terms of use.
 */

#ifndef NPGE_ALGO_LUA_HPP_
#define NPGE_ALGO_LUA_HPP_

#include <map>

#include "lua.hpp"
#include "luabind-format-signature.hpp"
#include <luabind/luabind.hpp>

#include "global.hpp"
#include "AnyAs.hpp"

namespace npge {

typedef std::map<std::string, AnyAs> MapAny;

}

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

template <>
struct default_converter<npge::MapAny> :
        native_converter_base<npge::MapAny> {
    static int compute_score(lua_State* L, int index);
    npge::MapAny from(lua_State* L, int index);
    void to(lua_State* L, const npge::MapAny& v);
};

template <>
struct default_converter<const npge::MapAny&> :
        default_converter<npge::MapAny> {
};

}

/** Initialize algo/ members in Lua */
extern "C" int init_algo_lua(lua_State* L);

#endif

