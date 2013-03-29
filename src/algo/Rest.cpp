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

Rest::Rest(const BlockSetPtr& source):
    skip_rest_(false) {
    set_other(source);
}

void Rest::add_options_impl(po::options_description& desc) const {
    add_unique_options(desc)
    ("skip-rest", po::value<bool>()->default_value(skip_rest_),
     "do not add unique fragments to block set")
   ;
}

void Rest::apply_options_impl(const po::variables_map& vm) {
    if (vm.count("skip-rest")) {
        skip_rest_ = vm["skip-rest"].as<bool>();
    }
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
            BOOST_ASSERT(!(*new_f < **prev));
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

bool Rest::run_impl() const {
    if (skip_rest_) {
        return false;
    }
    int blocks_before = block_set()->size();
    int seqs_before = block_set()->seqs().size();
    block_set()->add_sequences(other()->seqs());
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
                    BOOST_ASSERT(!(*f < *fr));
                    f = fr;
                }
                try_new_block(new_blocks, *f, -1, &prev);
                while (Fragment* fr = f->neighbor(1)) {
                    BOOST_ASSERT(!(*fr < *f));
                    f = fr;
                    try_new_block(new_blocks, *f, -1, &prev);
                }
                try_new_block(new_blocks, *f, 1, &prev);
            } else {
                BOOST_ASSERT(f->next() || f->prev());
            }
        }
    }
    size_t other_after = other()->size();
    BOOST_ASSERT(other_before == other_after);
    BOOST_FOREACH (Block* block, new_blocks) {
        block_set()->insert(block);
    }
    BOOST_FOREACH (SequencePtr seq, other()->seqs()) {
        if (used.find(seq.get()) == used.end()) {
            used.insert(seq.get());
            Block* block = new Block;
            block->set_name(seq->name());
            block->insert(new Fragment(seq, 0, seq->size() - 1));
            block_set()->insert(block);
        }
    }
    return block_set()->size() != blocks_before ||
           block_set()->seqs().size() != seqs_before;
}

}

