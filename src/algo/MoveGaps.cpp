/*
 * NPG-explorer, Nucleotide PanGenome explorer
 * Copyright (C) 2012-2015 Boris Nagaev
 *
 * See the LICENSE file for terms of use.
 */

#include <utility>
#include <boost/foreach.hpp>

#include "MoveGaps.hpp"
#include "Block.hpp"
#include "Fragment.hpp"
#include "AlignmentRow.hpp"
#include "RowStorage.hpp"
#include "throw_assert.hpp"
#include "cast.hpp"

namespace npge {

MoveGaps::MoveGaps() {
    add_row_storage_options(this);
    add_gopt("max-tail", "Max length of tail", "MAX_TAIL");
    add_gopt("max-tail-to-gap",
             "Max tail length to gap length ratio",
             "MAX_TAIL_TO_GAP");
    declare_bs("target", "Target blockset");
}

bool MoveGaps::move_gaps(Block* block) const {
    TimeIncrementer ti(this);
    int max_tail = opt_value("max-tail").as<int>();
    Decimal max_tail_to_gap = opt_value("max-tail-to-gap").as<Decimal>();
    int length = block->alignment_length();
    bool result = false;
    BOOST_FOREACH (Fragment* f, *block) {
        AlignmentRow* row = f->row();
        ASSERT_MSG(row, ("No alignment row is set, fragment " +
                         f->id()).c_str());
        ASSERT_MSG(row->length() == length,
                   ("Length of row of fragment " + f->id() +
                    " (" + TO_S(f->row()->length()) + ") "
                    "differs from block alignment length"
                    " (" + TO_S(length) + ")").c_str());
        typedef std::pair<int, int> TailGap;
        TailGap moves[3]; // index is ori + 1
        for (int ori = -1; ori <= 1; ori += 2) {
            int begin = (ori == 1) ? 0 : length - 1;
            int tail = 0;
            int i;
            for (i = 0; i < max_tail + 1; i++) {
                int al_pos = begin + i * ori;
                if (row->map_to_fragment(al_pos) != -1) {
                    tail += 1;
                } else {
                    break;
                }
            }
            if (0 < tail && tail <= max_tail) {
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
                    if (Decimal(tail) / gap <= max_tail_to_gap) {
                        moves[ori + 1] = std::make_pair(tail, gap);
                    }
                }
            }
        }
        if (moves[-1 + 1].first != 0 || moves[1 + 1].first != 0) {
            result = true;
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
            AlignmentRow* row = create_row(this);
            row->grow(data);
            f->set_row(row);
        }
    }
    return result;
}

void MoveGaps::process_block_impl(Block* block, ThreadData*) const {
    move_gaps(block);
}

const char* MoveGaps::name_impl() const {
    return "Move terminal letters inside";
}

}

