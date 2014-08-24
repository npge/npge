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
    static int compute_score(lua_State* L, int index) {
        return lua_type(L, index) == LUA_TTABLE ? 0 : -1;
    }

    // http://www.lua.org/manual/5.1/manual.html#lua_next
    npge::Strings from(lua_State* L, int index) {
        npge::Strings result;
        lua_pushnil(L); // first key
        while (lua_next(L, index) != 0) {
            // uses 'key' (at index -2)
            // and 'value' (at index -1)
            lua_getglobal(L, "tostring");
            lua_pushvalue(L, -2); // value
            lua_call(L, 1, 1);
            result.push_back(lua_tostring(L, -1));
            lua_pop(L, 2); // remove two 'value's
            // keep 'key' for next iteration
        }
        return result;
    }

    void to(lua_State* L, const npge::Strings& strings) {
        lua_createtable(L, strings.size(), 0);
        for (int i = 0; i < strings.size(); i++) {
            const std::string& s = strings[i];
            lua_pushstring(L, s.c_str());
            lua_rawseti(L, -2, i + 1);
        }
    }
};

template <>
struct default_converter<const npge::Strings&> :
        default_converter<npge::Strings> {
};

}

/** Initialize util/ members in Lua */
extern "C" int init_util_lua(lua_State* L);

#endif

