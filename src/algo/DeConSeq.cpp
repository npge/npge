/*
 * bloomrepeats, Find genomic repeats, using Bloom filter based prefiltration
 * Copyright (C) 2012 Boris Nagaev
 *
 * See the LICENSE file for terms of use.
 */

#include <boost/foreach.hpp>

#include "DeConSeq.hpp"
#include "Sequence.hpp"
#include "BlockSet.hpp"
#include "Block.hpp"
#include "Fragment.hpp"
#include "convert_position.hpp"
#include "throw_assert.hpp"

namespace bloomrepeats {

DeConSeq::DeConSeq(const BlockSetPtr& source) {
    set_other(source);
}

bool DeConSeq::run_impl() const {
    BOOST_FOREACH (const Block* block, *other()) {
        Block* new_block = new Block;
        new_block->set_name(block->name());
        block_set()->insert(new_block);
        BOOST_FOREACH (const Fragment* fragment, *block) {
            Sequence* seq = fragment->seq();
            BOOST_ASSERT(seq);
            const Block* seq_block = seq->block();
            BOOST_ASSERT_MSG(seq, "Sequence must be a consensus of a block");
            int seq_block_length = seq_block->alignment_length();
            BOOST_FOREACH (const Fragment* seq_f, *seq_block) {
                size_t fr_begin = fragment_pos(seq_f, fragment->begin_pos(),
                                               seq_block_length);
                size_t fr_last = fragment_pos(seq_f, fragment->last_pos(),
                                              seq_block_length);
                size_t seq_begin = frag_to_seq(seq_f, fr_begin);
                size_t seq_last = frag_to_seq(seq_f, fr_last);
                Fragment* new_fragment = new Fragment(seq_f->seq());
                new_fragment->set_begin_last(seq_begin, seq_last);
                new_block->insert(new_fragment);
            }
        }
    }
    return !other()->empty();
}

const char* DeConSeq::name_impl() const {
    return "Build blocks from blocks, set to sequences of blocks of other";
}

}

