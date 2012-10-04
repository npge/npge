/*
 * bloomrepeats, Find genomic repeats, using Bloom filter based prefiltration
 * Copyright (C) 2012 Boris Nagaev
 *
 * See the LICENSE file for terms of use.
 */

#include <boost/assert.hpp>
#include <boost/foreach.hpp>

#include "Joiner.hpp"
#include "Fragment.hpp"
#include "Block.hpp"

namespace bloomrepeats {

Joiner::Joiner(int max_dist, float ratio_to_fragment, float gap_ratio):
    max_dist_(max_dist),
    ratio_to_fragment_(ratio_to_fragment),
    gap_ratio_(gap_ratio)
{ }

bool Joiner::can_join_fragments(Fragment* f1, Fragment* f2) {
    BOOST_ASSERT(Fragment::can_join(f1, f2));
    int dist = f1->dist_to(*f2);
    int min_length = std::min(f1->length(), f2->length());
    BOOST_ASSERT(min_length > 0);
    float ratio = float(dist) / float(min_length);
    return (max_dist_ == -1 || dist <= max_dist_) &&
           (ratio_to_fragment_ < 0 || ratio <= ratio_to_fragment_);
}

bool Joiner::can_join_blocks(Block* b1, Block* b2) {
    BOOST_ASSERT(Block::can_join(b1, b2));
    BOOST_ASSERT(!b1->empty() && !b2->empty());
    Fragment* neighbor_1 = b1->front()->logical_neighbor(1);
    int ori = (neighbor_1 && neighbor_1->block() == b2) ? 1 : -1;
    BOOST_ASSERT(b1->front()->logical_neighbor(ori)->block() == b2);
    int min_gap = -1, max_gap = -1;
    BOOST_FOREACH (Fragment* f1, *b1) {
        Fragment* f2 = f1->logical_neighbor(ori);
        if (!can_join_fragments(f1, f2)) {
            return false;
        }
        int dist = f1->dist_to(*f2);
        min_gap = (min_gap == -1 || dist < min_gap) ? dist : min_gap;
        max_gap = (max_gap == -1 || dist > max_gap) ? dist : max_gap;
    }
    if (gap_ratio_ >= 0 && float(max_gap) / float(min_gap) > gap_ratio_) {
        return false;
    }
    return true;
}

}

