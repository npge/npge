/*
 * NPG-explorer, Nucleotide PanGenome explorer
 * Copyright (C) 2014 Boris Nagaev
 *
 * See the LICENSE file for terms of use.
 */

#include <iostream>
#include <string>
#include <boost/filesystem.hpp>

#include "process.hpp"
#include "meta_pipe.hpp"
#include "Processor.hpp"
#include "Meta.hpp"
#include "tss_meta.hpp"
#include "read_file.hpp"
#include "read_config.hpp"
#include "string_arguments.hpp"

using namespace npge;

int main(int argc, char** argv) {
    std::string app = argv[0];
    std::string script;
    if (argc >= 2 && argv[1][0] != '-') {
#if BOOST_FILESYSTEM_VERSION == 3
        app = boost::filesystem::path(argv[1]).filename().string();
#else
        app = boost::filesystem::path(argv[1]).filename();
#endif
        using namespace boost::filesystem;
        if (!exists(argv[1])) {
            std::cerr << "No such file: " << argv[1];
            std::cerr << std::endl;
            return 255;
        }
        script = read_file(argv[1]);
    }
    bool has_script = !script.empty();
    StringToArgv args(has_script ? app.c_str() : argv[0]);
    for (int i = has_script ? 2 : 1; i < argc; i++) {
        args.add_argument(argv[i]);
    }
    bool is_help = args.has_argument("-h") ||
                   args.has_argument("--help");
    if (argc == 1 || (!has_script && is_help)) {
        Processor p;
        print_help(":cout", &p, app, "");
        std::cout << std::endl;
        std::cout << "Pass script or '-i' as first argument";
        std::cout << std::endl;
        return 0;
    }
    bool interactive = args.has_argument("-i");
    args.remove_argument("-i");
    Meta& meta = *tss_meta();
    std::string c = args.get_argument("-c");
    if (!c.empty()) {
        meta.set_opt("LOCAL_CONF", c);
    }
    read_config(&meta);
    if (args.has_argument("-g")) {
        std::string g = args.get_argument("-g");
        if (g.empty()) {
            g = ":cout";
        }
        print_config(g, &meta);
        return 0;
    }
    int result = 0;
    bool debug = args.has_argument("--debug");
    if (has_script) {
        int r = execute_script(script, ":cerr",
                               args.argc(), args.argv(),
                               &meta, debug);
        if (r) {
            result = r;
        }
    }
    if (!interactive) {
        // non-interactive
        return result;
    } else {
        int r = interactive_loop(":cin", ":cout",
                                 args.argc(), args.argv(), &meta);
        return r ? : result;
    }
}

