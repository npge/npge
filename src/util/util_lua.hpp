/*
 * NPG-explorer, Nucleotide PanGenome explorer
 * Copyright (C) 2014 Boris Nagaev
 *
 * See the LICENSE file for terms of use.
 */

#ifndef NPGE_UTIL_LUA_HPP_
#define NPGE_UTIL_LUA_HPP_

#include "lua.hpp"
#include <luabind/luabind.hpp>

#include "global.hpp"

namespace luabind {

template <>
struct default_converter<npge::Strings> :
        native_converter_base<npge::Strings> {
    static int compute_score(lua_State* L, int index);
    npge::Strings from(lua_State* L, int index);
    void to(lua_State* L, const npge::Strings& strings);
};

template <>
struct default_converter<const npge::Strings&> :
        default_converter<npge::Strings> {
};

}

/** Initialize util/ members in Lua */
extern "C" int init_util_lua(lua_State* L);

#endif

