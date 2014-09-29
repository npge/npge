/*
 * NPG-explorer, Nucleotide PanGenome explorer
 * Copyright (C) 2014 Boris Nagaev
 *
 * See the LICENSE file for terms of use.
 */

#include <algorithm>
#include <boost/foreach.hpp>

#include "FragmentsExtender.hpp"
#include "MetaAligner.hpp"
#include "Sequence.hpp"
#include "AlignmentRow.hpp"
#include "Fragment.hpp"
#include "Block.hpp"
#include "complement.hpp"
#include "Decimal.hpp"
#include "throw_assert.hpp"

namespace npge {

FragmentsExtender::FragmentsExtender() {
    aligner_ = new MetaAligner;
    aligner_->set_parent(this);
    add_gopt("extend-length", "Length of extended part",
             "MIN_LENGTH");
    add_opt("extend-length-portion",
            "Length of extended part portion in source",
            D(0.0));
    declare_bs("target", "Target blockset");
}

static int max_right_shift(Fragment* f) {
    if (f->ori() == 1) {
        return int(f->seq()->size()) - 1 - int(f->max_pos());
    } else {
        return int(f->min_pos());
    }
}

static int max_right_shift(Block* block) {
    if (block->empty()) {
        return 0;
    }
    int result = max_right_shift(block->front());
    BOOST_FOREACH (Fragment* f, *block) {
        result = std::min(result, max_right_shift(f));
    }
    return result;
}

static void shift_end(Fragment* f, int shift) {
    if (f->ori() == 1) {
        f->set_max_pos(f->max_pos() + shift);
    } else {
        f->set_min_pos(f->min_pos() - shift);
    }
}

typedef std::map<Fragment*, std::string> F2S;

void extend_right(Block* block, F2S& right,
                  int extend_length, MetaAligner* aligner) {
    int right_length = max_right_shift(block);
    right_length = std::min(right_length, extend_length);
    if (right_length == 0) {
        return;
    }
    ASSERT_GT(right_length, 0);
    Fragments ff((block->begin()), block->end());
    Strings seqs;
    BOOST_FOREACH (Fragment* f, ff) {
        int start = f->length();
        int stop = start + right_length - 1;
        seqs.push_back(f->substr(start, stop));
        shift_end(f, right_length);
    }
    aligner->align_seqs(seqs);
    ASSERT_EQ(seqs.size(), ff.size());
    for (int i = 0; i < ff.size(); i++) {
        Fragment* f = ff[i];
        right[f].swap(seqs[i]);
    }
}

void FragmentsExtender::extend(Block* block) const {
    if (block->size() < 2 || !block->front()->row()) {
        // small or no alignment
        return;
    }
    F2S central;
    BOOST_FOREACH (Fragment* f, *block) {
        central[f] = f->str();
    }
    int extend_length = opt_value("extend-length").as<int>();
    Decimal portion;
    portion = opt_value("extend-length-portion").as<Decimal>();
    int length = block->alignment_length();
    int portion_length = (portion * length).to_i();
    extend_length = std::max(extend_length, portion_length);
    F2S right;
    extend_right(block, right, extend_length, aligner_);
    block->inverse(/* inverse_row */ false);
    F2S left;
    extend_right(block, left, extend_length, aligner_);
    block->inverse(/* inverse_row */ false);
    BOOST_FOREACH (Fragment* f, *block) {
        std::string& l = left[f];
        complement(l);
        const std::string& c = central[f];
        const std::string& r = right[f];
        AlignmentRow* row = f->row();
        ASSERT_TRUE(row);
        AlignmentRow* new_row = AlignmentRow::new_row(row->type());
        f->set_row(new_row);
        new_row->grow(l + c + r);
    }
}

void FragmentsExtender::process_block_impl(Block* block,
        ThreadData*) const {
    extend(block);
}

const char* FragmentsExtender::name_impl() const {
    return "Move block's boundaries and align only new parts";
}

}

