/*
 * NPG-explorer, Nucleotide PanGenome explorer
 * Copyright (C) 2012-2016 Boris Nagaev
 *
 * See the LICENSE file for terms of use.
 */

#include <iostream>
#include <vector>
#include <boost/foreach.hpp>
#include <boost/algorithm/string/predicate.hpp>

#include <lua.hpp>

#include "Meta.hpp"
#include "name_to_stream.hpp"

using namespace npge;

int main(int argc, char** argv) {
    if (argc < 2) {
        std::cerr << "Provide directory with tests" << "\n";
        return 255;
    }
    Meta meta;
    lua_State* L = meta.L();
    std::string test_dir = system_complete(argv[1]);
    if (!file_exists(test_dir)) {
        std::cerr << "Not found: " << test_dir << std::endl;
        return 255;
    }
    if (!is_dir(test_dir)) {
        std::cerr << "Not directory: " << test_dir << std::endl;
        return 255;
    }
    int all_scripts = 0;
    int ok_scripts = 0;
    Strings tests = dir_children(test_dir);
    BOOST_FOREACH (std::string test, tests) {
        using namespace boost::algorithm;
        if (!ends_with(test, ".lua")) {
            continue;
        }
        all_scripts += 1;
        int status = luaL_dofile(L, test.c_str());
        if (status) {
            std::cerr << "Error in file " << test << ":\n";
            std::cerr << lua_tostring(L, -1) << "\n";
            std::cerr << "\n";
        } else {
            ok_scripts += 1;
        }
    }
    std::cerr << "Total: " << all_scripts << " tests.\n";
    std::cerr << "Passed: " << ok_scripts << " tests.\n";
    std::cerr << std::endl;
    int failed_scripts = all_scripts - ok_scripts;
    std::cerr << "*** ";
    if (failed_scripts) {
        std::cerr << failed_scripts;
        std::cerr << " failed scripts detected";
        std::cerr << std::endl;
        return 255;
    } else {
        std::cerr << "No errors detected";
        std::cerr << std::endl;
        return 0;
    }
}

