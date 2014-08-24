/*
 * NPG-explorer, Nucleotide PanGenome explorer
 * Copyright (C) 2014 Boris Nagaev
 *
 * See the LICENSE file for terms of use.
 */

#include <iostream>
#include <lua.hpp>

#include "model_lua.hpp"
#include "util_lua.hpp"
#include "is_wine.hpp"
#include "terminal.lua"

#ifdef LUAPROMPT
extern "C" {
#include <term.h>
#include "luaprompt/prompt.h"
}
#endif

using namespace npge;

int main(int argc, char** argv) {
    lua_State* L = luaL_newstate();
    luaL_openlibs(L);
    init_model_lua(L);
    init_util_lua(L);
    if (argc >= 2) {
        // FIXME
        luaL_dofile(L, argv[1]);
        std::cerr << lua_tostring(L, -1) << "\n";
    }
    bool luaprompt = false;
#ifdef LUAPROMPT
    if (!is_wine()) {
        luaprompt = true;
        setupterm((char*)0, 1, (int*)0);
        tigetnum("colors");
        if (tigetnum("colors") <= 2) {
            luap_setcolor(L, 0);
        }
        luap_setname(L, "npge");
        luap_sethistory(L, "~/.npge_history");
        luap_enter(L);
    }
#endif
    if (!luaprompt) {
        // fallback terminal
        std::cerr << luaL_dostring(L, terminal_lua) << "\n";
        std::cerr << lua_tostring(L, -1) << "\n";
    }
    lua_close(L);
    std::cerr << "bye\n";
}

