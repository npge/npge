/*
 * bloomrepeats, Find genomic repeats, using Bloom filter based prefiltration
 * Copyright (C) 2012 Boris Nagaev
 *
 * See the LICENSE file for terms of use.
 */

#include <vector>
#include <boost/foreach.hpp>
#include <boost/cast.hpp>

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
    set_block_set_name("other");
}

struct CSData : public ThreadData {
    std::vector<SequencePtr> seqs_;
};

ThreadData* ConSeq::before_thread_impl() const {
    return new CSData;
}

void ConSeq::process_block_impl(Block* b, ThreadData* d) const {
    SequencePtr seq = create_sequence(this);
    seq->set_block(b);
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

