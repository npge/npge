/*
 * bloomrepeats, Find genomic repeats, using Bloom filter based prefiltration
 * Copyright (C) 2012 Boris Nagaev
 *
 * See the LICENSE file for terms of use.
 */

#include <algorithm>
#include <boost/foreach.hpp>

#include "OriByMajority.hpp"
#include "Block.hpp"
#include "Fragment.hpp"
#include "throw_assert.hpp"

namespace bloomrepeats {

static struct FragmentCompare3 {
    bool operator()(const Fragment* f1, const Fragment* f2) const {
        return *f1 < *f2;
    }
} fc3;

bool OriByMajority::apply_to_block_impl(Block* block) const {
    int ori_1 = 0;
    int sum = 0;
    BOOST_FOREACH (Fragment* f, *block) {
        sum += f->length();
        if (f->ori() == 1) {
            ori_1 += f->length();
        }
    }
    if (sum == 0) {
        return false;
    } else if (ori_1 * 2 < sum) {
        block->inverse();
        return true;
    } else if (ori_1 * 2 == sum) {
        Block::const_iterator it = std::min_element(block->begin(),
                                   block->end(), fc3);
        BOOST_ASSERT(it != block->end());
        Fragment* f = *it;
        if (f->ori() == -1) {
            block->inverse();
            BOOST_FOREACH (Fragment* f, *block) {
                f->set_row(0);
            }
            return true;
        }
    }
    return false;
}

const char* OriByMajority::name_impl() const {
    return "Set ori so that most nucleotides have ori=1";
}

}

