/*
 * bloomrepeats, Find genomic repeats, using Bloom filter based prefiltration
 * Copyright (C) 2012 Boris Nagaev
 *
 * See the LICENSE file for terms of use.
 */

#include <iostream>
#include <vector>
#include <boost/filesystem/path.hpp>
#include <boost/filesystem/operations.hpp>
#include <boost/foreach.hpp>

#include "process.hpp"
#include "meta_pipe.hpp"
#include "Meta.hpp"
#include "AddBlocks.hpp"
#include "name_to_stream.hpp"
#include "read_file.hpp"
#include "string_arguments.hpp"
#include "block_hash.hpp"

using namespace bloomrepeats;
namespace fs = boost::filesystem;

uint32_t hash_block_sets(const std::string& filename) {
    AddBlocks ab;
    ab.set_opt_value("in-blocks", filename);
    ab.run();
    uint32_t result = 0;
    Strings block_sets;
    ab.get_block_sets(block_sets);
    BOOST_FOREACH (const std::string& bs_name, block_sets) {
        result ^= blockset_hash(*ab.get_bs(bs_name));
    }
    return result;
}

bool run_test(const std::string& in_filename,
              const std::string& script_filename,
              const std::string& out_filename,
              const std::string& tmp_filename) {
    std::string script = read_file(script_filename);
    StringToArgv args;
    args.add_argument("--in-blocks");
    args.add_argument(in_filename);
    args.add_argument("--out-file");
    args.add_argument(tmp_filename);
    Meta meta;
    SharedProcessor p(parse_script(script, &meta));
    if (!p) {
        std::cerr << "Error: no processor found in "
                  << script_filename << std::endl;
        return false;
    }
    int r = process(args.argc(), args.argv(), p.get(), p->name());
    if (r != 0) {
        std::cerr << "Error executing " << script_filename << std::endl;
        std::cerr << " input file " << in_filename << std::endl;
        std::cerr << "Error code " << r << std::endl;
        return false;
    }
    uint32_t expected_hash = hash_block_sets(out_filename);
    uint32_t actual_hash = hash_block_sets(tmp_filename);
    if (expected_hash != actual_hash) {
        std::string out_expected = read_file(out_filename);
        std::string out_actual = read_file(tmp_filename);
        std::cerr << "Wrong output of " << script_filename << std::endl;
        std::cerr << "Input file: " << in_filename << std::endl;
        std::cerr << "Expected output:" << std::endl;
        std::cerr << out_expected << std::endl;
        std::cerr << "Actual output:" << std::endl;
        std::cerr << out_actual << std::endl;
        std::cerr << std::endl;
        return false;
    } else if (expected_hash == 0) {
        std::cerr << "Warning! Empty block set." << std::endl;
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
    fs::path test_dir = fs::system_complete(fs::path(argv[1]));
    if (!fs::exists(test_dir)) {
        std::cerr << "Not found: " << test_dir << std::endl;
        return 255;
    }
    if (!fs::is_directory(test_dir)) {
        std::cerr << "Not directory: " << test_dir << std::endl;
        return 255;
    }
    int all_scripts = 0, ok_scripts = 0;
    int all_tests = 0, ok_tests = 0;
    fs::directory_iterator dir(test_dir), end;
    BOOST_FOREACH (const fs::path& child_dir, std::make_pair(dir, end)) {
        if (fs::is_directory(child_dir)) {
            all_scripts += 1;
            bool script_ok = true;
            std::string script_filename = (child_dir / "script.br").string();
            fs::directory_iterator dir2(child_dir);
            BOOST_FOREACH (const fs::path& subtest, std::make_pair(dir2, end)) {
                if (fs::is_directory(subtest)) {
                    all_tests += 1;
                    std::string in_filename = (subtest / "in.fasta").string();
                    std::string out_filename = (subtest / "out.fasta").string();
                    std::string tmp_filename = ":test";
                    set_sstream(tmp_filename);
                    if (run_test(in_filename, script_filename,
                                 out_filename, tmp_filename)) {
                        ok_tests += 1;
                    } else {
                        script_ok = false;
                    }
                    remove_stream(tmp_filename);
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
    if (ok_tests != all_tests) {
        return 255;
    } else {
        return 0;
    }
}

