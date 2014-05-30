/*
 * NPG-explorer, Nucleotide PanGenome explorer
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
#include "Connector.hpp"
#include "throw_assert.hpp"

namespace npge {

Rest::Rest(const BlockSetPtr& source) {
    set_other(source);
    add_opt("skip-rest", "do not add unique fragments to block set",
            false);
    declare_bs("target", "Where created blocks are added");
    declare_bs("other", "Input blocks");
}

static void try_new_block(std::vector<Block*>& new_blocks,
                          const Fragment& f, int ori,
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
            ASSERT_FALSE(*new_f < **prev);
            Fragment::connect(*prev, new_f);
        }
        *prev = new_f;
        Block* block = new Block();
        block->insert(new_f);
        new_blocks.push_back(block);
    } else {
        delete new_f;
    }
}

void Rest::run_impl() const {
    if (opt_value("skip-rest").as<bool>()) {
        return;
    }
    Connector().apply(other());
    BlockSet& self = *block_set();
    self.add_sequences(other()->seqs());
    std::set<Sequence*> used;
    size_t other_before = other()->size();
    std::vector<Block*> new_blocks;
    BOOST_FOREACH (Block* block, *other()) {
        BOOST_FOREACH (Fragment* f, *block) {
            Sequence* seq = f->seq();
            if (used.find(seq) == used.end()) {
                used.insert(seq);
                Fragment* prev = 0;
                while (Fragment* fr = f->neighbor(-1)) {
                    ASSERT_FALSE(*f < *fr);
                    f = fr;
                }
                try_new_block(new_blocks, *f, -1, &prev);
                while (Fragment* fr = f->neighbor(1)) {
                    ASSERT_FALSE(*fr < *f);
                    f = fr;
                    try_new_block(new_blocks, *f, -1, &prev);
                }
                try_new_block(new_blocks, *f, 1, &prev);
            } else {
                ASSERT_TRUE(f->next() || f->prev());
            }
        }
    }
    size_t other_after = other()->size();
    ASSERT_EQ(other_before, other_after);
    BOOST_FOREACH (Block* block, new_blocks) {
        self.insert(block);
    }
    BOOST_FOREACH (SequencePtr seq, other()->seqs()) {
        if (used.find(seq.get()) == used.end()) {
            used.insert(seq.get());
            Block* block = new Block;
            block->set_name(seq->name());
            block->insert(new Fragment(seq, 0, seq->size() - 1));
            self.insert(block);
        }
    }
}

const char* Rest::name_impl() const {
    return "Add to target blocks of nucleotides, "
           "not included to other";
}

}

