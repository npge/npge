/*
 * NPG-explorer, Nucleotide PanGenome explorer
 * Copyright (C) 2014 Boris Nagaev
 *
 * See the LICENSE file for terms of use.
 */

#include <cstring>
#include <exception>
#include <iostream>
#include <boost/foreach.hpp>

#include "po.hpp"

namespace npge {

typedef boost::shared_ptr<po::option_description> OptPtr;

void add_new_options(const po::options_description& source,
                     po::options_description& dest,
                     const po::options_description* check) {
    BOOST_FOREACH (OptPtr opt, source.options()) {
        if (!dest.find_nothrow(opt->long_name(), /* approx */ false)) {
            if (!check || !check->find_nothrow(opt->long_name(), false)) {
                dest.add(opt);
            }
        }
    }
}

AddUniqueOptions::AddUniqueOptions(po::options_description& desc):
    po::options_description_easy_init(this),
    desc_(desc) {
}

AddUniqueOptions::~AddUniqueOptions() {
    add_new_options(*this, desc_);
}

AddUniqueOptions add_unique_options(po::options_description& desc) {
    return AddUniqueOptions(desc);
}

void add_general_options(po::options_description& desc) {
    add_unique_options(desc)
    ("help,h", "produce help message")
    ("debug", "do not catch errors")
    ("tree", "show processors tree")
    ("i", "interactive mode (tool npge)") // FIXME --i
    ("c", "path to local config file (option LOCAL_CONF)")
   ;
}

int read_options(int argc, char** argv, po::variables_map& vm,
                 const po::options_description& desc,
                 const po::positional_options_description& pod) {
    try {
        po::store(po::command_line_parser(argc, argv).
                  options(desc).positional(pod).run(), vm);
    } catch (std::exception& e) {
        std::cerr << argv[0] << ": error while parsing options: "
                  << std::endl << "  " << e.what() << std::endl;
        return 255;
    } catch (...) {
        std::cerr << argv[0] << ": error while parsing options: "
                  << std::endl << "  " << "Unknown error" << std::endl;
        return 255;
    }
    return 0;
}

bool has_arg(int argc, char** argv, const char* opt) {
    for (int i = 1; i < argc; i++) {
        if (strcmp(argv[i], opt) == 0) {
            return true;
        }
    }
    return false;
}

}

