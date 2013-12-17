/*
 * bloomrepeats, Find genomic repeats, using Bloom filter based prefiltration
 * Copyright (C) 2013 Boris Nagaev
 *
 * See the LICENSE file for terms of use.
 */

#include "hit.hpp"
#include "Block.hpp"
#include "Fragment.hpp"

namespace bloomrepeats {

bool is_internal_hit(const S2F& s2f, const Block* hit) {
    BOOST_FOREACH (Fragment* fragment, *hit) {
        std::vector<Fragment*> overlap_fragments;
        s2f.find_overlap_fragments(overlap_fragments, fragment);
        if (overlap_fragments.size() != 1) {
            return false; // multiple overlaps
        }
        Fragment* overlap_fragment = overlap_fragments[0];
        if (!fragment->is_subfragment_of(*overlap_fragment)) {
            return false; // hits is not inside block's fragment
        }
        Block* overlap_block = overlap_fragment->block();
        if (overlap_block->size() >= hit->size()) {
            return false; // block's size is >= than hit's size
        }
    }
    return true;
}

}

