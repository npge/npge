/*
 * bloomrepeats, Find genomic repeats, using Bloom filter based prefiltration
 * Copyright (C) 2014 Boris Nagaev
 *
 * See the LICENSE file for terms of use.
 */

#include <boost/foreach.hpp>

#include "MergeUnique.hpp"
#include "Connector.hpp"
#include "BlockSet.hpp"
#include "Block.hpp"
#include "Fragment.hpp"
#include "global.hpp"
#include "throw_assert.hpp"

namespace bloomrepeats {

MergeUnique::MergeUnique() {
    declare_bs("target", "Target blockset");
}

static void inspect_neighbours(Block* b, BlockSet& bs, int ori) {
    BOOST_ASSERT(b->size() >= 2);
    typedef std::pair<Block*, int> BlockOri;
    typedef std::map<BlockOri, Fragments> UniqueOf;
    UniqueOf unique_of;
    BOOST_FOREACH (Fragment* f, *b) {
        Fragment* n = f->logical_neighbor(ori);
        if (n && n->block() && n->block()->size() == 1) {
            Fragment* in_2 = (n == f->next()) ?
                             n->next() : n->prev();
            BOOST_ASSERT(in_2 != f);
            if (in_2) {
                Block* in_2_b = in_2->block();
                BOOST_ASSERT(in_2_b);
                BOOST_ASSERT(in_2_b->size() >= 2);
                int o = f->ori() * in_2->ori();
                unique_of[BlockOri(in_2_b, o)].push_back(n);
            }
        }
    }
    BOOST_FOREACH (const UniqueOf::value_type& u, unique_of) {
        const Fragments& ff = u.second;
        if (ff.size() >= 2) {
            Block* new_block = new Block;
            bs.insert(new_block);
            BOOST_FOREACH (Fragment* n, ff) {
                Block* n_b = n->block();
                new_block->insert(n);
                BOOST_ASSERT(n_b->size() == 1);
                n_b->detach(n);
                BOOST_ASSERT(n_b->size() == 0);
                bs.erase(n_b);
                Fragment* f = (n->prev()->block() == b) ?
                              n->prev() : n->next();
                BOOST_ASSERT(f->block() == b);
                n->set_ori(f->ori());
            }
        }
    }
}

void MergeUnique::run_impl() const {
    Connector c;
    c.apply(block_set());
    BlockSet& bs = *block_set();
    Blocks blocks;
    BOOST_FOREACH (Block* block, bs) {
        if (block->size() >= 2) {
            blocks.push_back(block);
        }
    }
    BOOST_FOREACH (Block* b, blocks) {
        for (int ori = -1; ori <= 1; ori += 2) {
            inspect_neighbours(b, bs, ori);
        }
    }
}

const char* MergeUnique::name_impl() const {
    return "Merge unique fragments with both common neighbours";
}

}

