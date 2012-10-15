/*
 * bloomrepeats, Find genomic repeats, using Bloom filter based prefiltration
 * Copyright (C) 2012 Boris Nagaev
 *
 * See the LICENSE file for terms of use.
 */

#include <iostream>
#include <fstream>
#include <boost/foreach.hpp>
#include <boost/algorithm/string/replace.hpp>

#include "Output.hpp"
#include "Exception.hpp"
#include "BlockSet.hpp"
#include "Block.hpp"

namespace bloomrepeats {

void Output::add_options_impl(po::options_description& desc) const {
    add_unique_options(desc)
    ("out-file,o", po::value<std::string>()->default_value(file()),
     "output file with all blocks")
    ("out-mask", po::value<std::string>()->default_value(mask()),
     "mask of output files (${block} is replaced with block name)");
}

void Output::apply_options_impl(const po::variables_map& vm) {
    set_file(vm["out-file"].as<std::string>());
    set_mask(vm["out-mask"].as<std::string>());
    if (!file().empty() && !mask().empty()) {
        throw Exception("both 'out-file' and 'out-mask' were specified");
    }
    if (!mask().empty() && mask().find("${block}") == std::string::npos) {
        throw Exception("'out-mask' must contain '${block}'");
    }
}

bool Output::run_impl() const {
    if (!mask().empty()) {
        BOOST_FOREACH (Block* b, *block_set()) {
            using namespace boost::algorithm;
            std::string path = replace_all_copy(mask(), "${block}", b->name());
            std::ofstream o(path.c_str());
            o << *b << std::endl;
        }
    }
    if (!file().empty()) {
        std::ofstream o(file().c_str());
        o << *block_set() << std::endl;
    }
    if (file().empty() && mask().empty()) {
        std::cout << *block_set() << std::endl;
    }
    return false;
}

const char* Output::name_impl() const {
    return "Output block set";
}

}

