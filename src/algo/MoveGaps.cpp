/*
 * bloomrepeats, Find genomic repeats, using Bloom filter based prefiltration
 * Copyright (C) 2012 Boris Nagaev
 *
 * See the LICENSE file for terms of use.
 */

#include <utility>
#include <boost/foreach.hpp>
#include <boost/lexical_cast.hpp>

#include "MoveGaps.hpp"
#include "Block.hpp"
#include "Fragment.hpp"
#include "AlignmentRow.hpp"
#include "throw_assert.hpp"

namespace bloomrepeats {

MoveGaps::MoveGaps(int max_tail, float max_tail_to_gap):
    max_tail_(max_tail), max_tail_to_gap_(max_tail_to_gap)
{ }

void MoveGaps::add_options_impl(po::options_description& desc) const {
    RowStorage::add_options_impl(desc);
    add_unique_options(desc)
    ("max-tail", po::value<int>()->default_value(max_tail()),
     "Max length of tail")
    ("max-tail-to-gap", po::value<float>()->default_value(max_tail_to_gap()),
     "Max tail length to gap length ratio")
   ;
}

void MoveGaps::apply_options_impl(const po::variables_map& vm) {
    RowStorage::apply_options_impl(vm);
    if (vm.count("max-tail")) {
        set_max_tail(vm["max-tail"].as<int>());
    }
    if (vm.count("max-tail-to-gap")) {
        set_max_tail_to_gap(vm["max-tail-to-gap"].as<float>());
    }
}

bool MoveGaps::apply_to_block_impl(Block* block) const {
    int length = block->alignment_length();
    BOOST_FOREACH (Fragment* f, *block) {
        AlignmentRow* row = f->row();
        BOOST_ASSERT_MSG(row, ("No alignment row is set, fragment " +
                               f->id()).c_str());
        BOOST_ASSERT_MSG(row->length() == length,
                         ("Length of row of fragment " + f->id() + " (" +
                          boost::lexical_cast<std::string>(f->row()->length()) +
                          ") differs from block alignment length (" +
                          boost::lexical_cast<std::string>(length)).c_str());
        typedef std::pair<int, int> TailGap;
        TailGap moves[3]; // index is ori + 1
        for (int ori = -1; ori <= 1; ori += 2) {
            int begin = (ori == 1) ? 0 : length - 1;
            int tail = 0;
            int i;
            for (i = 0; i < max_tail() + 1; i++) {
                int al_pos = begin + i * ori;
                if (row->map_to_fragment(al_pos) != -1) {
                    tail += 1;
                } else {
                    break;
                }
            }
            if (0 < tail && tail <= max_tail()) {
                int gap = 0;
                int max_pos = length / 2;
                for (; i < max_pos; i++) {
                    int al_pos = begin + i * ori;
                    if (row->map_to_fragment(al_pos) == -1) {
                        gap += 1;
                    } else {
                        break;
                    }
                }
                if (i < max_pos && gap != 0) {
                    if (float(tail) / float(gap) <= max_tail_to_gap()) {
                        moves[ori + 1] = std::make_pair(tail, gap);
                    }
                }
            }
        }
        if (moves[-1 + 1].first != 0 || moves[1 + 1].first != 0) {
            std::stringstream ss;
            f->print_contents(ss, '-', /* line */ 0);
            std::string data = ss.str();
            for (int ori = -1; ori <= 1; ori += 2) {
                int begin = (ori == 1) ? 0 : length - 1;
                int tail = moves[ori + 1].first;
                int gap = moves[ori + 1].second;
                if (tail) {
                    for (int i = tail - 1; i >= 0; i--) {
                        int to = gap + i;
                        int from = i;
                        data[begin + to * ori] = data[begin + from * ori];
                    }
                    for (int i = 0; i < gap; i++) {
                        data[begin + i * ori] = '-';
                    }
                }
            }
            AlignmentRow* row = AlignmentRow::new_row(row_type());
            row->grow(data);
            f->set_row(row);
        }
    }
}

const char* MoveGaps::name_impl() const {
    return "Move terminal letters inside";
}

}

