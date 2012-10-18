/*
 * bloomrepeats, Find genomic repeats, using Bloom filter based prefiltration
 * Copyright (C) 2012 Boris Nagaev
 *
 * See the LICENSE file for terms of use.
 */

#include <fstream>
#include <boost/foreach.hpp>

#include "AddBlocks.hpp"
#include "BlockSet.hpp"
#include "Alignment.hpp"

namespace bloomrepeats {

AddBlocks::AddBlocks()
{ }

AddBlocks::AddBlocks(Alignment* alignment):
    alignment_(alignment)
{ }

void AddBlocks::add_options_impl(po::options_description& desc) const {
    add_unique_options(desc)
    ("in-blocks", po::value<Files>()->required(),
     "input fasta file(s) with blocks")
   ;
}

void AddBlocks::apply_options_impl(const po::variables_map& vm) {
    set_files(vm["in-blocks"].as<Files>());
}

bool AddBlocks::run_impl() const {
    int size_before = block_set()->size();
    if (alignment() && alignment()->block_set() != block_set()) {
        alignment()->set_block_set(block_set());
    }
    BOOST_FOREACH (std::string file_name, files()) {
        std::ifstream input_file(file_name.c_str());
        if (alignment()) {
            input_file >> *alignment();
        } else {
            input_file >> *block_set();
        }
    }
    return block_set()->size() > size_before;
}

const char* AddBlocks::name_impl() const {
    return "Input block set";
}

}

