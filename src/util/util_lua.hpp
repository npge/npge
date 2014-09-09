/*
 * NPG-explorer, Nucleotide PanGenome explorer
 * Copyright (C) 2014 Boris Nagaev
 *
 * See the LICENSE file for terms of use.
 */

#ifndef NPGE_UTIL_LUA_HPP_
#define NPGE_UTIL_LUA_HPP_

#include <limits.h>
#include <stdint.h>

#include "lua.hpp"
#include "luabind-format-signature.hpp"
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

template <>
struct default_converter<npge::AnyAs> :
        native_converter_base<npge::AnyAs> {
    static int compute_score(lua_State* L, int index);
    npge::AnyAs from(lua_State* L, int index);
    void to(lua_State* L, const npge::AnyAs& strings);
};

template <>
struct default_converter<const npge::AnyAs&> :
        default_converter<npge::AnyAs> {
};

#ifdef LUABIND_INT64_MISSING

template <>
struct default_converter<int64_t> :
        native_converter_base<int64_t> {
    static int compute_score(lua_State* L, int index);
    int64_t from(lua_State* L, int index);
    void to(lua_State* L, const int64_t& value);
};

template <>
struct default_converter<const int64_t&> :
        default_converter<int64_t> {
};

#endif

}

namespace npge {

/** Set value of global var "arg" */
void set_arg(lua_State* L, const Strings& arg);

}

/** Initialize util/ members in Lua */
extern "C" int init_util_lua(lua_State* L);

#endif

