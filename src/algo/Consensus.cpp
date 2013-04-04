/*
 * bloomrepeats, Find genomic repeats, using Bloom filter based prefiltration
 * Copyright (C) 2012 Boris Nagaev
 *
 * See the LICENSE file for terms of use.
 */

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
    if (vm.count("out-consensus")) {
        set_output_file(vm["out-consensus"].as<std::string>());
    }
}

bool Consensus::run_impl() const {
    BOOST_FOREACH (Block* b, *block_set()) {
        output() << ">" << b->name();
        output() << std::endl;
        b->consensus(output());
        output() << std::endl;
    }
    return false;
}

const char* Consensus::name_impl() const {
    return "Consensus writer";
}

}

