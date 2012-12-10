/*
 * bloomrepeats, Find genomic repeats, using Bloom filter based prefiltration
 * Copyright (C) 2012 Boris Nagaev
 *
 * See the LICENSE file for terms of use.
 */

#include <fstream>
#include <boost/foreach.hpp>

#include "Consensus.hpp"
#include "Block.hpp"
#include "BlockSet.hpp"

namespace bloomrepeats {

void Consensus::add_options_impl(po::options_description& desc) const {
    add_unique_options(desc)
    ("out-consensus", po::value<std::string>()->required(),
     "Output file with consensuses")
   ;
}

void Consensus::apply_options_impl(const po::variables_map& vm) {
    set_file(vm["out-consensus"].as<std::string>());
}

bool Consensus::run_impl() const {
    std::ofstream out(file().c_str());
    BOOST_FOREACH (Block* b, *block_set()) {
        out << ">" << b->name();
        out << std::endl;
        b->consensus(out);
        out << std::endl;
    }
}

const char* Consensus::name_impl() const {
    return "Consensus writer";
}

}

