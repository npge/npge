/*
 * bloomrepeats, Find genomic repeats, using Bloom filter based prefiltration
 * Copyright (C) 2012 Boris Nagaev
 *
 * See the LICENSE file for terms of use.
 */

#include <cstdio>
#include <iostream>
#include <boost/filesystem/path.hpp>
#include <boost/filesystem/operations.hpp>
#include <boost/foreach.hpp>

#include "process.hpp"
#include "meta_pipe.hpp"
#include "Meta.hpp"
#include "temp_file.hpp"
#include "read_file.hpp"
#include "string_arguments.hpp"

using namespace bloomrepeats;
namespace fs = boost::filesystem;

bool run_test(const std::string& in_filename,
              const std::string& script_filename,
              const std::string& out_filename,
              const std::string& tmp_filename) {
    std::string script = read_file(script_filename);
    std::string out_reference = read_file(out_filename);
    StringToArgv args;
    args.add_argument("--in-blocks");
    args.add_argument(in_filename);
    args.add_argument("--out-file");
    args.add_argument(tmp_filename);
    Meta meta;
    ProcessorPtr p = parse_script(script, &meta);
    int r = process(args.argc(), args.argv(), p, p->name());
    if (r != 0) {
        std::cerr << "Error executing " << script_filename << std::endl;
        std::cerr << "Error code " << r << std::endl;
        return false;
    }
    std::string out_actual = read_file(tmp_filename);
    if (out_actual != out_reference) {
        std::cerr << "Wrong output of " << script_filename << std::endl;
        std::cerr << "Expected output:" << std::endl;
        std::cerr << out_reference << std::endl;
        std::cerr << "Actual output:" << std::endl;
        std::cerr << out_actual << std::endl;
        std::cerr << std::endl;
        return false;
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
    std::string tmp_filename = temp_file();
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
                    if (run_test(in_filename, script_filename,
                                 out_filename, tmp_filename)) {
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
    std::remove(tmp_filename.c_str());
    std::cerr << "Total: " << all_scripts << " test scripts, ";
    std::cerr << all_tests << " tests." << std::endl;
    std::cerr << "Passed: " << ok_scripts << " test scripts, ";
    std::cerr << ok_tests << " tests." << std::endl;
}

