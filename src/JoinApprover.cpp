/*
 * bloomrepeats, Find genomic repeats, using Bloom filter based prefiltration
 * Copyright (C) 2012 Boris Nagaev
 *
 * See the LICENSE file for terms of use.
 */

#include <boost/assert.hpp>
#include <boost/foreach.hpp>

#include "JoinApprover.hpp"
#include "Fragment.hpp"
#include "Block.hpp"

namespace bloomrepeats {

JoinApprover::JoinApprover(int max_dist):
    max_dist_(max_dist)
{ }

bool JoinApprover::can_join_fragments(FragmentPtr f1, FragmentPtr f2) {
    BOOST_ASSERT(Fragment::can_join(f1, f2));
    return max_dist_ == -1 || int(f1->dist_to(*f2)) <= max_dist_;
}

bool JoinApprover::can_join_blocks(BlockPtr b1, BlockPtr b2) {
    BOOST_ASSERT(Block::can_join(b1, b2));
    BOOST_ASSERT(!b1->empty() && !b2->empty());
    FragmentPtr neighbor_1 = b1->front()->logical_neighbor(1);
    int ori = (neighbor_1 && neighbor_1->block() == b2) ? 1 : -1;
    BOOST_ASSERT(b1->front()->logical_neighbor(ori)->block() == b2);
    BOOST_FOREACH (FragmentPtr f1, *b1) {
        FragmentPtr f2 = f1->logical_neighbor(ori);
        if (!can_join_fragments(f1, f2)) {
            return false;
        }
    }
    // TODO other checks
    return true;
}

}

