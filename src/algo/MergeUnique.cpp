/*
 * NPG-explorer, Nucleotide PanGenome explorer
 * Copyright (C) 2012-2016 Boris Nagaev
 *
 * See the LICENSE file for terms of use.
 */

#include <set>
#include <algorithm>
#include <boost/foreach.hpp>

#include "MergeUnique.hpp"
#include "FragmentCollection.hpp"
#include "BlockSet.hpp"
#include "Block.hpp"
#include "Fragment.hpp"
#include "block_hash.hpp"
#include "Meta.hpp"
#include "global.hpp"
#include "throw_assert.hpp"

namespace npge {

MergeUnique::MergeUnique() {
    add_opt("both-neighbours", "Require both neighbours "
            "of an unique fragment to be from common blocks "
            "(otherwise only one common neighbour is "
            "sufficient)", true);
    add_opt("merge-long", "Merge long fragments "
            "(>= MIN_LENGTH) as well", false);
    declare_bs("target", "Target blockset");
}

typedef std::set<Fragment*> FragmentsSet;

static void merge_fragments(
    const FragmentsSet& ff,
    const VectorFc& fc,
    Block* b,
    BlockSet& bs,
    int ori) {
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
        Fragment* f = (n_prev && n_prev->block() == b)
                      ? n_prev : fc.next(n);
        ASSERT_TRUE(f);
        ASSERT_EQ(f->block(), b);
        n->set_ori(f->ori());
    }
    new_block->set_name("m" + block_id(new_block));
}

static bool isJoinableFragment(const Fragment* n,
                               int min_length) {
    if (!n) {
        return false;
    }
    const Block* block = n->block();
    if (!block) {
        return false;
    }
    if (block->size() != 1) {
        return false;
    }
    if (block->alignment_length() >= min_length) {
        return false;
    }
    return true;
}

// merge unique fragments surrounded by same blocks
static void inspect_neighbours2(const VectorFc& fc,
                                Block* b, BlockSet& bs,
                                int ori, int min_length) {
    ASSERT_GTE(b->size(), 2);
    typedef std::pair<Block*, int> BlockOri;
    typedef std::map<BlockOri, Fragments> UniqueOf;
    UniqueOf unique_of;
    BOOST_FOREACH (Fragment* f, *b) {
        Fragment* n = fc.logical_neighbor(f, ori);
        if (isJoinableFragment(n, min_length)) {
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
        FragmentsSet ff;
        BOOST_FOREACH (Fragment* f, ff0) {
            ASSERT_TRUE(f->block());
            if (isJoinableFragment(f, min_length)) {
                ff.insert(f);
            }
        }
        if (ff.size() >= 2) {
            merge_fragments(ff, fc, b, bs, ori);
        }
    }
}

// merge unique neighbours of a block
static void inspect_neighbours1(const VectorFc& fc,
                                Block* b, BlockSet& bs,
                                int ori, int min_length) {
    ASSERT_GTE(b->size(), 2);
    FragmentsSet unique;
    BOOST_FOREACH (Fragment* f, *b) {
        Fragment* n = fc.logical_neighbor(f, ori);
        if (isJoinableFragment(n, min_length)) {
            unique.insert(n);
        }
    }
    if (unique.size() >= 2) {
        merge_fragments(unique, fc, b, bs, ori);
    }
}

struct BlockSizeCmpRev {
    bool operator()(const Block* b1, const Block* b2) const {
        return b2->size() < b1->size();
    }
};

void MergeUnique::run_impl() const {
    bool both = opt_value("both-neighbours").as<bool>();
    bool merge_long = opt_value("merge-long").as<bool>();
    int min_length = meta()->get_opt("MIN_LENGTH").as<int>();
    if (merge_long) {
        min_length = npge::MAX_POS;
    }
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
    std::sort(blocks.begin(), blocks.end(), BlockSizeCmpRev());
    BOOST_FOREACH (Block* b, blocks) {
        for (int ori = -1; ori <= 1; ori += 2) {
            if (both) {
                inspect_neighbours2(fc, b, bs, ori, min_length);
            } else {
                inspect_neighbours1(fc, b, bs, ori, min_length);
            }
        }
    }
}

const char* MergeUnique::name_impl() const {
    return "Merge unique fragments with "
           "common neighbours";
}

}

