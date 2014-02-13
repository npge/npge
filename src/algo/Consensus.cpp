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

Consensus::Consensus():
    file_writer_(this, "out-consensus", "Output file with consensus", true)
{ }

bool Consensus::run_impl() const {
    std::ostream& out = file_writer_.output();
    BOOST_FOREACH (Block* b, *block_set()) {
        out << ">" << b->name();
        out << std::endl;
        b->consensus(out);
        out << std::endl;
    }
    return false;
}

const char* Consensus::name_impl() const {
    return "Consensus writer";
}

}

