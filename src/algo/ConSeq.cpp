/*
 * bloomrepeats, Find genomic repeats, using Bloom filter based prefiltration
 * Copyright (C) 2012 Boris Nagaev
 *
 * See the LICENSE file for terms of use.
 */

#include <boost/foreach.hpp>

#include "ConSeq.hpp"
#include "SeqStorage.hpp"
#include "Sequence.hpp"
#include "BlockSet.hpp"

namespace bloomrepeats {

ConSeq::ConSeq(const BlockSetPtr& source) {
    set_other(source);
    add_seq_storage_options(this);
}

void ConSeq::add_options_impl(po::options_description& desc) const
{ }

void ConSeq::apply_options_impl(const po::variables_map& vm)
{ }

bool ConSeq::run_impl() const {
    bool result = false;
    BOOST_FOREACH (const Block* block, *other()) {
        SequencePtr seq = create_sequence(this);
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

