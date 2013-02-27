/*
 * bloomrepeats, Find genomic repeats, using Bloom filter based prefiltration
 * Copyright (C) 2012 Boris Nagaev
 *
 * See the LICENSE file for terms of use.
 */

#include <set>
#include <boost/foreach.hpp>

#include "Rest.hpp"
#include "Block.hpp"
#include "Fragment.hpp"
#include "Sequence.hpp"
#include "BlockSet.hpp"
#include "throw_assert.hpp"

namespace bloomrepeats {

Rest::Rest(const BlockSetPtr& source) {
    set_other(source);
}

static void try_new_block(BlockSet& set, const Fragment& f, int ori,
                          Fragment** prev) {
    Fragment* n = f.neighbor(ori);
    Fragment* new_f = new Fragment(f.seq());
    if (ori == -1) {
        new_f->set_min_pos(n ? n->max_pos() + 1 : 0);
        new_f->set_max_pos(f.min_pos() - 1);
    } else {
        new_f->set_min_pos(f.max_pos() + 1);
        new_f->set_max_pos(n ? n->min_pos() - 1 : f.seq()->size() - 1);
    }
    if (new_f->valid()) {
        if (*prev) {
            BOOST_ASSERT(!(*new_f < **prev));
            Fragment::connect(*prev, new_f);
        }
        *prev = new_f;
        Block* block = new Block();
        block->insert(new_f);
        set.insert(block);
    } else {
        delete new_f;
    }
}

bool Rest::run_impl() const {
    int blocks_before = block_set()->size();
    int seqs_before = block_set()->seqs().size();
    block_set()->add_sequences(other()->seqs());
    std::set<Sequence*> used;
    BOOST_FOREACH (Block* block, *other()) {
        BOOST_FOREACH (Fragment* f, *block) {
            Sequence* seq = f->seq();
            if (used.find(seq) == used.end()) {
                used.insert(seq);
                Fragment* prev = 0;
                while (Fragment* fr = f->neighbor(-1)) {
                    f = fr;
                }
                try_new_block(*block_set(), *f, -1, &prev);
                while (Fragment* fr = f->neighbor(1)) {
                    f = fr;
                    try_new_block(*block_set(), *f, -1, &prev);
                }
                try_new_block(*block_set(), *f, 1, &prev);
            }
        }
    }
    BOOST_FOREACH (SequencePtr seq, other()->seqs()) {
        if (used.find(seq.get()) == used.end()) {
            used.insert(seq.get());
            Block* block = new Block;
            block->insert(new Fragment(seq, 0, seq->size() - 1));
            block_set()->insert(block);
        }
    }
    return block_set()->size() != blocks_before ||
           block_set()->seqs().size() != seqs_before;
}

}

