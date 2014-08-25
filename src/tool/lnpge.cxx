/*
 * NPG-explorer, Nucleotide PanGenome explorer
 * Copyright (C) 2014 Boris Nagaev
 *
 * See the LICENSE file for terms of use.
 */

#include <iostream>
#include <lua.hpp>

#include "Meta.hpp"
#include "tss_meta.hpp"
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
    Meta& meta = *tss_meta();
    lua_State* L = meta.L();
    luaL_openlibs(L);
    if (argc >= 2) {
        // FIXME
        luaL_dofile(L, argv[1]);
        std::cerr << lua_tostring(L, -1) << "\n";
    }
    bool luaprompt = false;
#ifdef LUAPROMPT
    if (!is_wine()) {
        luaprompt = true;
#if defined(_WIN32) || defined(__WIN32__)
        luap_setcolor(L, 0);
#else
        setupterm((char*)0, 1, (int*)0);
        if (tigetnum("colors") <= 2) {
            luap_setcolor(L, 0);
        }
#endif
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
    std::cerr << "bye\n";
}

