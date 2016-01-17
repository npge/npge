/*
 * NPG-explorer, Nucleotide PanGenome explorer
 * Copyright (C) 2012-2016 Boris Nagaev
 *
 * See the LICENSE file for terms of use.
 */

#include <iostream>
#include <vector>
#include <boost/foreach.hpp>

#include "luabind-error.hpp"
#include "luabind-format-signature.hpp"
#include <lua.hpp>
#include <luabind/luabind.hpp>

#include "process.hpp"
#include "util_lua.hpp"
#include "Meta.hpp"
#include "Read.hpp"
#include "name_to_stream.hpp"
#include "read_file.hpp"
#include "string_arguments.hpp"
#include "block_hash.hpp"

using namespace npge;

hash_t hash_block_sets(const std::string& filename) {
    Read ab;
    ab.set_opt_value("in-blocks", filename);
    ab.run();
    hash_t result = 0;
    Strings block_sets;
    ab.get_block_sets(block_sets);
    BOOST_FOREACH (const std::string& bs_name, block_sets) {
        result ^= blockset_hash(*ab.get_bs(bs_name));
    }
    return result;
}

struct RemoveStream {
    std::string file_;

    RemoveStream(const std::string& file):
        file_(file) {
    }

    ~RemoveStream() {
        remove_stream(file_);
    }
};

bool run_test(const std::string& in_filename,
              const std::string& script_filename,
              const std::string& out_filename,
              Meta& meta) {
    lua_State* L = meta.L();
    std::string tmp_filename = ":test";
    set_sstream(tmp_filename);
    RemoveStream rm(tmp_filename);
    StringToArgv args;
    args.add_argument(script_filename);
    args.add_argument("--in-blocks");
    args.add_argument(in_filename);
    args.add_argument("--out-file");
    args.add_argument(tmp_filename);
    set_arg(L, args.to_strings());
    meta.reset_placeholder_processor();
    int r = luaL_dostring(L, "main()");
    if (r) {
        std::cerr << "Error executing " << script_filename;
        std::cerr << std::endl;
        std::cerr << " input file " << in_filename << std::endl;
        std::cerr << "Error code " << r << std::endl;
        std::cerr << lua_tostring(L, -1) << "\n";
        return false;
    }
    std::string out_actual = read_file(tmp_filename);
    set_sstream(tmp_filename, out_actual);
    hash_t expected_hash = hash_block_sets(out_filename);
    hash_t actual_hash = hash_block_sets(tmp_filename);
    if (expected_hash != actual_hash) {
        std::string out_expected = read_file(out_filename);
        std::cerr << "Wrong output of " << script_filename << std::endl;
        std::cerr << "Input file: " << in_filename << std::endl;
        std::cerr << "Expected output:" << std::endl;
        std::cerr << out_expected << std::endl;
        std::cerr << "Actual output:" << std::endl;
        std::cerr << out_actual << std::endl;
        std::cerr << std::endl;
        return false;
    } else if (expected_hash == 0) {
        std::cerr << "Warning! Empty blockset." << std::endl;
        std::cerr << "Output file: " << out_filename << std::endl;
        std::cerr << std::endl;
    }
    return true;
}

int main(int argc, char** argv) {
    if (argc < 2) {
        std::cerr << "Provide directory with tests" << std::endl;
        return 255;
    }
    Meta meta;
    std::string test_dir = system_complete(argv[1]);
    if (!file_exists(test_dir)) {
        std::cerr << "Not found: " << test_dir << std::endl;
        return 255;
    }
    if (!is_dir(test_dir)) {
        std::cerr << "Not directory: " << test_dir << std::endl;
        return 255;
    }
    int all_scripts = 0, ok_scripts = 0;
    int all_tests = 0, ok_tests = 0;
    Strings tests = dir_children(test_dir);
    BOOST_FOREACH (std::string child_dir, tests) {
        if (is_dir(child_dir)) {
            bool script_ok = true;
            std::string script_filename;
            script_filename = cat_paths(child_dir,
                                        "script.npge");
            if (!file_exists(script_filename)) {
                continue;
            }
            all_scripts += 1;
            Strings subtests = dir_children(child_dir);
            BOOST_FOREACH (std::string subtest, subtests) {
                if (is_dir(subtest)) {
                    all_tests += 1;
                    std::string in_filename, out_filename;
                    in_filename = cat_paths(subtest,
                                            "in.fasta");
                    out_filename = cat_paths(subtest,
                                             "out.fasta");
                    if (run_test(in_filename, script_filename,
                                 out_filename, meta)) {
                        ok_tests += 1;
                    } else {
                        script_ok = false;
                    }
                }
            }
            if (script_ok) {
                ok_scripts += 1;
            }
        }
    }
    std::cerr << "Total: " << all_scripts << " test scripts, ";
    std::cerr << all_tests << " tests." << std::endl;
    std::cerr << "Passed: " << ok_scripts << " test scripts, ";
    std::cerr << ok_tests << " tests." << std::endl;
    std::cerr << std::endl;
    int failed_scripts = all_scripts - ok_scripts;
    int failed_tests = all_tests - ok_tests;
    std::cerr << "*** ";
    if (failed_scripts) {
        std::cerr << failed_tests << " failed tests in ";
        std::cerr << failed_scripts << " scripts detected";
        std::cerr << std::endl;
        return 255;
    } else {
        std::cerr << "No errors detected";
        std::cerr << std::endl;
        return 0;
    }
}

