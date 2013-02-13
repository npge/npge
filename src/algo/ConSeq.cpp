/*
 * bloomrepeats, Find genomic repeats, using Bloom filter based prefiltration
 * Copyright (C) 2012 Boris Nagaev
 *
 * See the LICENSE file for terms of use.
 */

#include <boost/foreach.hpp>

#include "ConSeq.hpp"
#include "Sequence.hpp"
#include "BlockSet.hpp"

namespace bloomrepeats {

ConSeq::ConSeq(const BlockSetPtr& source) {
    set_other(source);
}

void ConSeq::add_options_impl(po::options_description& desc) const {
    SeqStorage::add_options_impl(desc);
}

void ConSeq::apply_options_impl(const po::variables_map& vm) {
    SeqStorage::apply_options_impl(vm);
}

bool ConSeq::run_impl() const {
    bool result = false;
    BOOST_FOREACH (const Block* block, *other()) {
        SequencePtr seq = create_sequence();
        seq->set_block(block);
        block_set()->add_sequence(seq);
        result = true;
    }
    return result;
}

const char* ConSeq::name_impl() const {
    return "Blocks to consensus sequences";
}

}

