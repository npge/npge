/*
 * bloomrepeats, Find genomic repeats, using Bloom filter based prefiltration
 * Copyright (C) 2013 Boris Nagaev
 *
 * See the LICENSE file for terms of use.
 */

#include <boost/foreach.hpp>
#include <boost/tuple/tuple.hpp>
#include <boost/tuple/tuple_comparison.hpp>

#include "hit.hpp"
#include "Block.hpp"
#include "Fragment.hpp"
#include "convert_position.hpp"
#include "throw_assert.hpp"

namespace bloomrepeats {

bool is_internal_hit(const S2F& s2f, const Block* hit,
                     bool allow_no_overlaps) {
    BOOST_FOREACH (Fragment* fragment, *hit) {
        std::vector<Fragment*> overlap_fragments;
        s2f.find_overlap_fragments(overlap_fragments, fragment);
        if (overlap_fragments.empty() && allow_no_overlaps) {
            continue;
        }
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

class FragmentCompareSequenceFirst {
public:
    bool operator()(const Fragment* a, const Fragment* b) const {
        typedef boost::tuple<const Sequence*, const Fragment&> Tie;
        return Tie(a->seq(), *a) < Tie(b->seq(), *b);
    }
};

bool has_self_overlaps(Block* block) {
    if (block->empty()) {
        return false;
    }
    std::vector<Fragment*> fragments(block->begin(), block->end());
    std::sort(fragments.begin(), fragments.end(),
              FragmentCompareSequenceFirst());
    for (int i = 0; i < fragments.size() - 1; i++) {
        Fragment* current = fragments[i];
        Fragment* next = fragments[i + 1];
        if (current->common_positions(*next)) {
            return true;
        }
    }
    return false;
}

void fix_self_overlaps(Block* block) {
    if (!has_self_overlaps(block)) {
        return;
    }
    boost::shared_ptr<Block> copy((block->clone()));
    int block_length = copy->alignment_length();
    for (int length = block_length - 1; length >= 0; length--) {
        block->clear();
        BOOST_FOREACH (Fragment* f, *copy) {
            int fragment_last_pos = fragment_pos(f, length, block_length);
            size_t seq_last = frag_to_seq(f, fragment_last_pos);
            size_t seq_begin = f->begin_pos();
            if (seq_last != seq_begin) {
                Fragment* new_f = new Fragment(f->seq());
                new_f->set_begin_last(seq_begin, seq_last);
                block->insert(new_f);
            }
        }
        if (!has_self_overlaps(block)) {
            break;
        }
    }
    BOOST_ASSERT(!has_self_overlaps(block));
}

}

