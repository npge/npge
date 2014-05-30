/*
 * NPG-explorer, Nucleotide PanGenome explorer
 * Copyright (C) 2014 Boris Nagaev
 *
 * See the LICENSE file for terms of use.
 */

#include <boost/foreach.hpp>
#include <boost/scoped_ptr.hpp>
#include <boost/cast.hpp>

#include "DeConSeq.hpp"
#include "Sequence.hpp"
#include "BlockSet.hpp"
#include "Block.hpp"
#include "Fragment.hpp"
#include "AlignmentRow.hpp"
#include "convert_position.hpp"
#include "throw_assert.hpp"
#include "Exception.hpp"

namespace npge {

DeConSeq::DeConSeq(const BlockSetPtr& source) {
    set_other(source);
    declare_bs("target", "Where new blocks are added");
    declare_bs("other", "From where consensus blocks are taken");
    set_block_set_name("other");
}

static void deconseq_row(const Fragment* fragment, Fragment* f) {
    AlignmentRow* row = fragment->row();
    boost::scoped_ptr<AlignmentRow> tmp_row((f->detach_row()));
    ASSERT_TRUE(row);
    ASSERT_TRUE(tmp_row);
    AlignmentRow* new_row = AlignmentRow::new_row(row->type());
    f->set_row(new_row);
    int length = row->length();
    new_row->set_length(length);
    for (int i = 0; i < length; i++) {
        int tmp_pos = row->map_to_fragment(i);
        if (tmp_pos != -1) {
            int new_pos = tmp_row->map_to_fragment(tmp_pos);
            if (new_pos != -1) {
                new_row->bind(new_pos, i);
            }
        }
    }
}

Block* DeConSeq::deconseq_block(const Block* block) {
    Block* new_block = new Block;
    new_block->set_name(block->name());
    BOOST_FOREACH (const Fragment* fragment, *block) {
        Sequence* seq = fragment->seq();
        ASSERT_TRUE(seq);
        const Block* seq_block = seq->block();
        if (!seq_block) {
            throw Exception("Sequence " + seq->name() +
                            " is not bound to a block");
        }
        int start = fragment->begin_pos();
        int stop = fragment->last_pos();
        bool alignment = (fragment->row() != 0);
        boost::scoped_ptr<Block> temp_block((seq_block->slice(start,
                                             stop, alignment)));
        BOOST_FOREACH (Fragment* f, *temp_block) {
            new_block->insert(f); // fragment->block() == new_block
            if (f->row()) {
                ASSERT_TRUE(fragment->row());
                deconseq_row(fragment, f);
            }
        }
        // temp_block is deleted, but new fragments do not
    }
    return new_block;
}

struct DSData : public ThreadData {
    Blocks blocks_;
};

ThreadData* DeConSeq::before_thread_impl() const {
    return new DSData;
}

void DeConSeq::process_block_impl(Block* b,
                                  ThreadData* d) const {
    DSData* data = boost::polymorphic_cast<DSData*>(d);
    Block* new_block = deconseq_block(b);
    data->blocks_.push_back(new_block);
}

void DeConSeq::after_thread_impl(ThreadData* d) const {
    BlockSet& t = *block_set();
    DSData* data = boost::polymorphic_cast<DSData*>(d);
    BOOST_FOREACH (Block* b, data->blocks_) {
        t.insert(b);
    }
}

const char* DeConSeq::name_impl() const {
    return "Build blocks from blocks, set to sequences of blocks of other";
}

}

