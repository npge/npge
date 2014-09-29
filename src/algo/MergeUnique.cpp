/*
 * NPG-explorer, Nucleotide PanGenome explorer
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
#include "block_hash.hpp"
#include "global.hpp"
#include "throw_assert.hpp"

namespace npge {

MergeUnique::MergeUnique() {
    declare_bs("target", "Target blockset");
}

static void inspect_neighbours(Block* b, BlockSet& bs, int ori) {
    ASSERT_GTE(b->size(), 2);
    typedef std::pair<Block*, int> BlockOri;
    typedef std::map<BlockOri, Fragments> UniqueOf;
    UniqueOf unique_of;
    BOOST_FOREACH (Fragment* f, *b) {
        Fragment* n = f->logical_neighbor(ori);
        if (n && n->block() && n->block()->size() == 1) {
            Fragment* in_2 = (n == f->next()) ?
                             n->next() : n->prev();
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
                Fragment* f = (n->prev()->block() == b) ?
                              n->prev() : n->next();
                ASSERT_EQ(f->block(), b);
                n->set_ori(f->ori());
            }
            new_block->set_name("m" + block_id(new_block));
        }
    }
}

void MergeUnique::run_impl() const {
    Connector c;
    c.set_opt_value("connect-circular", true);
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
    c.set_opt_value("connect-circular", false);
    c.apply(block_set());
}

const char* MergeUnique::name_impl() const {
    return "Merge unique fragments with both common neighbours";
}

}

