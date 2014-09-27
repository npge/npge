/*
 * NPG-explorer, Nucleotide PanGenome explorer
 * Copyright (C) 2014 Boris Nagaev
 *
 * See the LICENSE file for terms of use.
 */

#include <iostream>

#include "luabind-format-signature.hpp"
#include <lua.hpp>
#include <luabind/luabind.hpp>
#include <luabind/tag_function.hpp>

#include "Meta.hpp"
#include "opts_lib.hpp"
#include "process.hpp"
#include "is_wine.hpp"
#include "name_to_stream.hpp"
#include "string_arguments.hpp"
#include "util_lua.hpp"
#include "model_lua.hpp"
#include "algo_lua.hpp"
#include "npge_debug.hpp"

#ifdef LUAPROMPT
extern "C" {
#include <term.h>
#include "luaprompt/prompt.h"
}
#endif

using namespace npge;

void npge_terminal(lua_State* L);

int main(int argc, char** argv) {
    StringToArgv args(argc, argv);
    std::string c = args.get_argument("-c");
    if (!c.empty()) {
        set_local_conf(c);
    }
    Meta meta;
    if (meta.get_opt("NPGE_DEBUG").as<bool>()) {
        set_npge_debug(true);
    }
    lua_State* L = meta.L();
    set_arg(L, args.to_strings());
    using namespace luabind;
    module(L) [
        def("terminal", tag_function<void()>(
                boost::bind(npge_terminal, L)))
    ];
    int status = luaL_dostring(L, "main()");
    if (status) {
        std::cerr << lua_tostring(L, -1) << "\n";
    }
    return status;
}

void npge_terminal(lua_State* L) {
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
        std::string hist = get_home_dir() + "/.npge_history";
        luap_sethistory(L, hist.c_str());
        luap_enter(L);
    }
#endif
    if (!luaprompt) {
        luaL_dostring(L, "simple_terminal()");
    }
}

