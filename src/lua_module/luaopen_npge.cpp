/*
 * NPG-explorer, Nucleotide PanGenome explorer
 * Copyright (C) 2014 Boris Nagaev
 *
 * See the LICENSE file for terms of use.
 */

#include "luaopen_npge.hpp"
#include "Meta.hpp"

extern "C" int luaopen_npge(lua_State* L) {
    npge::Meta* meta = new npge::Meta; // FIXME memory leak?
    meta->attach_to_lua(L);
    lua_getglobal(L, "meta");
    return 1; // return meta
}

