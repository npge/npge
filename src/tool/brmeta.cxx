/*
 * bloomrepeats, Find genomic repeats, using Bloom filter based prefiltration
 * Copyright (C) 2012 Boris Nagaev
 *
 * See the LICENSE file for terms of use.
 */

#include <iostream>
#include <string>
#include <boost/filesystem.hpp>

#include "process.hpp"
#include "meta_pipe.hpp"
#include "Meta.hpp"
#include "read_file.hpp"
#include "string_arguments.hpp"

using namespace bloomrepeats;

int main(int argc, char** argv) {
    if (argc < 2) {
        std::cerr << "Pass script as first argument" << std::endl;
        return 255;
    }
    std::string script = read_file(argv[1]);
#if BOOST_FILESYSTEM_VERSION == 3
    std::string app = boost::filesystem::path(argv[1]).filename().string();
#else
    std::string app = boost::filesystem::path(argv[1]).filename();
#endif
    StringToArgv args(app.c_str());
    for (int i = 2; i < argc; i++) {
        args.add_argument(argv[i]);
    }
    Meta meta;
    ProcessorPtr p = parse_script(script, &meta);
    if (!p) {
        std::cerr << "Can't find main processor in script" << std::endl;
        return 255;
    }
    return process(args.argc(), args.argv(), p, p->name());
}

