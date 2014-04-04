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
    declare_bs("other",
               "Blockset, from which blocks are taken");
    declare_bs("target",
               "Blockset, where consensus sequences are added");
}

void ConSeq::run_impl() const {
    BOOST_FOREACH (const Block* block, *other()) {
        SequencePtr seq = create_sequence(this);
        seq->set_block(block);
        block_set()->add_sequence(seq);
    }
}

const char* ConSeq::name_impl() const {
    return "Blocks to consensus sequences";
}

}

