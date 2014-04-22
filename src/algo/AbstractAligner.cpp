/*
 * bloomrepeats, Find genomic repeats, using Bloom filter based prefiltration
 * Copyright (C) 2014 Boris Nagaev
 *
 * See the LICENSE file for terms of use.
 */

#include <boost/foreach.hpp>

#include "AbstractAligner.hpp"
#include "AlignmentRow.hpp"
#include "Block.hpp"
#include "Fragment.hpp"
#include "throw_assert.hpp"
#include "to_s.hpp"

namespace bloomrepeats {

AbstractAligner::AbstractAligner() {
    declare_bs("target", "Target blockset");
}

struct BlockSquareLess {
    bool operator()(Block* a, Block* b) const {
        int a_length = a->front() ? a->front()->length() : 0;
        int b_length = b->front() ? b->front()->length() : 0;
        return b_length * b->size() < a_length * a->size();
    }
};

void AbstractAligner::change_blocks_impl(Blocks& blocks) const {
    std::sort(blocks.begin(), blocks.end(), BlockSquareLess());
}

void AbstractAligner::align_block(Block* block) const {
    if (!alignment_needed(block)) {
        return;
    }
    Fragments fragments((block->begin()), block->end());
    Strings rows;
    BOOST_FOREACH (Fragment* f, fragments) {
        rows.push_back(f->str(/* gap */ 0));
    }
    align_seqs(rows);
    ASSERT_EQ(rows.size(), fragments.size());
    for (int i = 0; i < fragments.size(); i++) {
        new CompactAlignmentRow(rows[i], fragments[i]);
    }
}

void AbstractAligner::align_seqs(Strings& seqs) const {
    int size_before = seqs.size();
    align_seqs_impl(seqs);
    int size_after = seqs.size();
    ASSERT_EQ(size_after, size_before);
    if (size_after > 0) {
        int length = seqs.front().length();
        BOOST_FOREACH (const std::string& seq, seqs) {
            ASSERT_EQ(seq.length(), length);
        }
    }
}

bool AbstractAligner::alignment_needed(Block* block) const {
    if (block->size() == 0) {
        return false;
    } else if (block->size() == 1) {
        Fragment* f = block->front();
        if (f->row() && f->row()->length() == f->length()) {
            return false;
        }
        AlignmentRow* row = AlignmentRow::new_row(COMPACT_ROW);
        int length = f->length();
        row->set_length(length);
        for (int i = 0; i < length; i++) {
            row->bind(i, i);
        }
        f->set_row(row);
        return false;
    }
    if (block->front()->row()) {
        int row_length = block->front()->row()->length();
        bool all_rows = true;
        BOOST_FOREACH (Fragment* f, *block) {
            if (!f->row() || f->row()->length() != row_length) {
                all_rows = false;
                break;
            }
        }
        if (all_rows) {
            // all fragments have rows and lengthes are equal
            return false;
        }
    }
    return true;
}

void AbstractAligner::process_block_impl(Block* block,
        ThreadData*) const {
    align_block(block);
}

const char* AbstractAligner::name_impl() const {
    return "Align blocks";
}

}

