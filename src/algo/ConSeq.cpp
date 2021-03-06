/*
 * NPG-explorer, Nucleotide PanGenome explorer
 * Copyright (C) 2012-2016 Boris Nagaev
 *
 * See the LICENSE file for terms of use.
 */

#include <vector>
#include <boost/make_shared.hpp>
#include <boost/foreach.hpp>
#include <boost/cast.hpp>

#include "ConSeq.hpp"
#include "SeqStorage.hpp"
#include "Sequence.hpp"
#include "BlockSet.hpp"
#include "Block.hpp"

namespace npge {

ConSeq::ConSeq(const BlockSetPtr& source) {
    set_other(source);
    add_seq_storage_options(this);
    declare_bs("other",
               "Blockset, from which blocks are taken");
    declare_bs("target",
               "Blockset, where consensus sequences are added");
    set_block_set_name("other");
}

struct CSData : public ThreadData {
    std::vector<SequencePtr> seqs_;
};

ThreadData* ConSeq::before_thread_impl() const {
    return new CSData;
}

void ConSeq::process_block_impl(Block* b, ThreadData* d) const {
    SequencePtr seq;
    if (b->size() == 1) {
        seq = boost::make_shared<FragmentSequence>(b->front());
        seq->set_block(b, /* set consensus */ false);
    } else if (b->size() >= 2) {
        seq = create_sequence(this);
        seq->set_block(b);
    }
    CSData* data = boost::polymorphic_cast<CSData*>(d);
    data->seqs_.push_back(seq);
}

void ConSeq::after_thread_impl(ThreadData* d) const {
    BlockSet& t = *block_set();
    CSData* data = boost::polymorphic_cast<CSData*>(d);
    BOOST_FOREACH (const SequencePtr& seq, data->seqs_) {
        t.add_sequence(seq);
    }
}

const char* ConSeq::name_impl() const {
    return "Blocks to consensus sequences";
}

}

