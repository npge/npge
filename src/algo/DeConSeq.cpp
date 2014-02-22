/*
 * bloomrepeats, Find genomic repeats, using Bloom filter based prefiltration
 * Copyright (C) 2012 Boris Nagaev
 *
 * See the LICENSE file for terms of use.
 */

#include <boost/foreach.hpp>
#include <boost/scoped_ptr.hpp>

#include "DeConSeq.hpp"
#include "Sequence.hpp"
#include "BlockSet.hpp"
#include "Block.hpp"
#include "Fragment.hpp"
#include "convert_position.hpp"
#include "throw_assert.hpp"
#include "Exception.hpp"

namespace bloomrepeats {

DeConSeq::DeConSeq(const BlockSetPtr& source) {
    set_other(source);
}

Block* DeConSeq::deconseq_block(const Block* block) {
    Block* new_block = new Block;
    new_block->set_name(block->name());
    BOOST_FOREACH (const Fragment* fragment, *block) {
        Sequence* seq = fragment->seq();
        BOOST_ASSERT(seq);
        const Block* seq_block = seq->block();
        if (!seq_block) {
            throw Exception("Sequence " + seq->name() +
                            " is not bound to a block");
        }
        int start = fragment->begin_pos();
        int stop = fragment->last_pos();
        boost::scoped_ptr<Block> temp_block((seq_block->slice(start, stop)));
        BOOST_FOREACH (Fragment* f, *temp_block) {
            new_block->insert(f); // fragment->block() == new_block
        }
        // temp_block is deleted, but new fragments do not
    }
    return new_block;
}

bool DeConSeq::run_impl() const {
    BOOST_FOREACH (const Block* block, *other()) {
        Block* new_block = deconseq_block(block);
        block_set()->insert(new_block);
    }
    return !other()->empty();
}

const char* DeConSeq::name_impl() const {
    return "Build blocks from blocks, set to sequences of blocks of other";
}

}

