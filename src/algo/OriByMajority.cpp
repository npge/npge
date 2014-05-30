/*
 * NPG-explorer, Nucleotide PanGenome explorer
 * Copyright (C) 2014 Boris Nagaev
 *
 * See the LICENSE file for terms of use.
 */

#include <algorithm>
#include <boost/foreach.hpp>
#include <boost/tuple/tuple.hpp>
#include <boost/tuple/tuple_comparison.hpp>

#include "OriByMajority.hpp"
#include "Block.hpp"
#include "Fragment.hpp"
#include "Sequence.hpp"
#include "throw_assert.hpp"

namespace npge {

OriByMajority::OriByMajority() {
    declare_bs("target", "Target blockset");
}

static struct FragmentCompare3 {
    bool operator()(const Fragment* f1, const Fragment* f2) const {
        typedef boost::tuple<size_t, size_t, const std::string&> Tie;
        ASSERT_TRUE(f1->seq());
        ASSERT_TRUE(f2->seq());
        return Tie(f1->min_pos(), f1->max_pos(), f1->seq()->name()) <
               Tie(f2->min_pos(), f2->max_pos(), f2->seq()->name());
    }
} fc3;

void OriByMajority::process_block_impl(Block* block, ThreadData*) const {
    bool result = false;
    int ori_1 = 0;
    int sum = 0;
    BOOST_FOREACH (Fragment* f, *block) {
        sum += f->length();
        if (f->ori() == 1) {
            ori_1 += f->length();
        }
    }
    if (sum == 0) {
        result = false;
    } else if (ori_1 * 2 < sum) {
        result = true;
    } else if (ori_1 * 2 == sum) {
        Block::const_iterator it = std::min_element(block->begin(),
                                   block->end(), fc3);
        ASSERT_TRUE(it != block->end());
        Fragment* f = *it;
        if (f->ori() == -1) {
            result = true;
        }
    }
    if (result == true) {
        block->inverse();
    }
}

const char* OriByMajority::name_impl() const {
    return "Set ori so that most nucleotides have ori=1";
}

}

