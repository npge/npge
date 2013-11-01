/*
 * bloomrepeats, Find genomic repeats, using Bloom filter based prefiltration
 * Copyright (C) 2012 Boris Nagaev
 *
 * See the LICENSE file for terms of use.
 */

#include <sstream>
#include <algorithm>
#include <boost/foreach.hpp>
#include <boost/tuple/tuple.hpp>
#include <boost/tuple/tuple_comparison.hpp>

#include "OriByMajority.hpp"
#include "Block.hpp"
#include "Fragment.hpp"
#include "Sequence.hpp"
#include "AlignmentRow.hpp"
#include "complement.hpp"
#include "throw_assert.hpp"

namespace bloomrepeats {

static struct FragmentCompare3 {
    bool operator()(const Fragment* f1, const Fragment* f2) const {
        typedef boost::tuple<size_t, size_t, const std::string&> Tie;
        BOOST_ASSERT(f1->seq()->name());
        BOOST_ASSERT(f2->seq()->name());
        return Tie(f1->min_pos(), f1->max_pos(), f1->seq()->name()) <
               Tie(f2->min_pos(), f2->max_pos(), f2->seq()->name());
    }
} fc3;

bool OriByMajority::apply_to_block_impl(Block* block) const {
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
        BOOST_ASSERT(it != block->end());
        Fragment* f = *it;
        if (f->ori() == -1) {
            result = true;
        }
    }
    if (result == true) {
        block->inverse();
        BOOST_FOREACH (Fragment* f, *block) {
            AlignmentRow* row = f->row();
            if (row) {
                std::stringstream ss;
                f->print_contents(ss, '-', /* line */ 0);
                std::string data = ss.str();
                complement(data);
                AlignmentRow* new_row = AlignmentRow::new_row(COMPACT_ROW);
                // TODO row type
                new_row->grow(data);
                f->set_row(new_row);
            }
        }
    }
    return result;
}

const char* OriByMajority::name_impl() const {
    return "Set ori so that most nucleotides have ori=1";
}

}

