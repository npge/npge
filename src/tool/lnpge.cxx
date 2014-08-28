/*
 * NPG-explorer, Nucleotide PanGenome explorer
 * Copyright (C) 2014 Boris Nagaev
 *
 * See the LICENSE file for terms of use.
 */

#include <iostream>
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

#ifdef LUAPROMPT
extern "C" {
#include <term.h>
#include "luaprompt/prompt.h"
}
#endif

using namespace npge;

void lnpge_terminal(lua_State* L);

int main(int argc, char** argv) {
    std::string app = argv[0];
    set_app_path(app);
    bool has_script = (argc >= 2 && argv[1][0] != '-');
    StringToArgv args(has_script ? app.c_str() : argv[0]);
    for (int i = has_script ? 2 : 1; i < argc; i++) {
        args.add_argument(argv[i]);
    }
    std::string c = args.get_argument("-c");
    if (!c.empty()) {
        set_local_conf(c);
    }
    Meta meta;
    if (args.has_argument("-g")) {
        std::string g = args.get_argument("-g");
        if (g.empty()) {
            g = ":cout";
        }
        print_config(g, &meta);
        return 0;
    }
    lua_State* L = meta.L();
    luabind::globals(L)["main_args"] = args.to_s();
    if (has_script) {
        int status = luaL_dofile(L, argv[1]);
        if (status) {
            std::cerr << lua_tostring(L, -1) << "\n";
        }
        if (!args.has_argument("-i")) {
            return status;
        }
    }
    using namespace luabind;
    module(L) [
        def("terminal", tag_function<void()>(
            boost::bind(lnpge_terminal, L)))
    ];
    int status = luaL_dostring(L, "terminal()");
    if (status) {
        std::cerr << lua_tostring(L, -1) << "\n";
    }
}

void lnpge_terminal(lua_State* L) {
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
    std::cerr << "bye\n";
}

