/*
 * NPG-explorer, Nucleotide PanGenome explorer
 * Copyright (C) 2014 Boris Nagaev
 *
 * See the LICENSE file for terms of use.
 */

#include <boost/foreach.hpp>

#include "MergeUnique.hpp"
#include "FragmentCollection.hpp"
#include "BlockSet.hpp"
#include "Block.hpp"
#include "Fragment.hpp"
#include "block_hash.hpp"
#include "global.hpp"
#include "throw_assert.hpp"

namespace npge {

MergeUnique::MergeUnique() {
    declare_bs("target", "Target blockset");
}

static void inspect_neighbours(const VectorFc& fc,
                               Block* b, BlockSet& bs,
                               int ori) {
    ASSERT_GTE(b->size(), 2);
    typedef std::pair<Block*, int> BlockOri;
    typedef std::map<BlockOri, Fragments> UniqueOf;
    UniqueOf unique_of;
    BOOST_FOREACH (Fragment* f, *b) {
        Fragment* n = fc.logical_neighbor(f, ori);
        if (n && n->block() && n->block()->size() == 1) {
            Fragment* in_2 = fc.another_neighbor(n, f);
            // in_2 can be == f,
            // if the sequence has only 2 fragments
            if (in_2 && in_2 != f) {
                Block* in_2_b = in_2->block();
                ASSERT_TRUE(in_2_b);
                if (in_2_b->size() >= 2) {
                    int o = f->ori() * in_2->ori();
                    unique_of[BlockOri(in_2_b, o)].push_back(n);
                }
            }
        }
    }
    BOOST_FOREACH (const UniqueOf::value_type& u, unique_of) {
        const Fragments& ff0 = u.second;
        // exclude used fragments
        typedef std::set<Fragment*> FragmentsSet;
        FragmentsSet ff;
        BOOST_FOREACH (Fragment* f, ff0) {
            ASSERT_TRUE(f->block());
            if (f->block()->size() == 1) {
                ff.insert(f);
            }
        }
        if (ff.size() >= 2) {
            Block* new_block = new Block;
            bs.insert(new_block);
            BOOST_FOREACH (Fragment* n, ff) {
                Block* n_b = n->block();
                new_block->insert(n);
                ASSERT_EQ(n_b->size(), 1);
                n_b->detach(n);
                ASSERT_EQ(n_b->size(), 0);
                bs.erase(n_b);
                Fragment* n_prev = fc.prev(n);
                Fragment* f = (n_prev->block() == b)
                              ? n_prev : fc.next(n);
                ASSERT_EQ(f->block(), b);
                n->set_ori(f->ori());
            }
            new_block->set_name("m" + block_id(new_block));
        }
    }
}

void MergeUnique::run_impl() const {
    BlockSet& bs = *block_set();
    VectorFc fc;
    fc.add_bs(bs);
    fc.prepare();
    Blocks blocks;
    BOOST_FOREACH (Block* block, bs) {
        if (block->size() >= 2) {
            blocks.push_back(block);
        }
    }
    BOOST_FOREACH (Block* b, blocks) {
        for (int ori = -1; ori <= 1; ori += 2) {
            inspect_neighbours(fc, b, bs, ori);
        }
    }
}

const char* MergeUnique::name_impl() const {
    return "Merge unique fragments with both "
           "common neighbours";
}

}

