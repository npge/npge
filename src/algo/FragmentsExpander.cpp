/*
 * bloomrepeats, Find genomic repeats, using Bloom filter based prefiltration
 * Copyright (C) 2012 Boris Nagaev
 *
 * See the LICENSE file for terms of use.
 */

#include <boost/foreach.hpp>

#include "FragmentsExpander.hpp"
#include "ExpanderBase.hpp"
#include "PairAligner.hpp"
#include "Fragment.hpp"
#include "Block.hpp"
#include "BlockSet.hpp"
#include "Exception.hpp"

namespace bloomrepeats {

static struct BlockGreater2 {
    bool operator()(const Block* b1, const Block* b2) const {
        return b1->size() > b2->size();
    }
} block_greater_2;

FragmentsExpander::FragmentsExpander(int batch, int ori,
                                     int max_overlap) {
    add_expander_options(this);
    set_opt_value("batch", batch);
    add_opt("ori", "direction of expansion", ori);
    add_opt("max-overlap",
            "max number of positions added after first overlap occur",
            max_overlap);
    add_opt_rule("max-overlap >= 0");
}

bool FragmentsExpander::change_blocks_impl(std::vector<Block*>& bs) const {
    std::sort(bs.begin(), bs.end(), block_greater_2);
    return false;
}

bool FragmentsExpander::process_block_impl(Block* block, ThreadData*) const {
    bool result = expand(block);
    if (result) {
        BOOST_FOREACH (Fragment* f, *block) {
            f->set_row(0);
        }
    }
    return result;
}

const char* FragmentsExpander::name_impl() const {
    return "Expand fragments";
}

bool FragmentsExpander::expand(Block* block) const {
    if (block->size() < 2) {
        return false;
    }
    int max_errors = opt_value("max-errors").as<int>();
    int gap_range = opt_value("gap-range").as<int>();
    int gap_penalty = opt_value("gap-penalty").as<int>();
    int ori = opt_value("ori").as<int>();
    PairAligner aligner_copy(max_errors, gap_range, gap_penalty);
    bool result = false;
    if (ori == 1 || ori == 0) {
        result |= expand_end(block, aligner_copy);
    }
    if (ori == -1 || ori == 0) {
        block->inverse();
        result |= expand_end(block, aligner_copy);
        block->inverse();
    }
    return result;
}

bool FragmentsExpander::expand_end(Block* block, PairAligner& a) const {
    int batch = opt_value("batch").as<int>();
    int max_overlap = opt_value("max-overlap").as<int>();
    if (block->size() <= 1) {
        return false;
    }
    std::vector<int> main_end(block->size() - 1), o_end(block->size() - 1);
    Fragment* main_f = block->front();
    bool result = false;
    while (true) {
        int max_shift = block->max_shift_end(max_overlap);
        if (max_shift <= 0) {
            break;
        }
        result = true;
        int shift = std::min(batch, max_shift);
        std::string main_str = main_f->substr(-1, main_f->length() - 1 + shift);
        a.set_first(main_str.c_str(), main_str.size());
        int i = 0;
        BOOST_FOREACH (Fragment* o_f, *block) {
            if (o_f != main_f) {
                std::string o_str = o_f->substr(-1, o_f->length() - 1 + shift);
                a.set_second(o_str.c_str(), o_str.size());
                a.align(main_end[i], o_end[i]);
                i++;
            }
        }
        int min_end = *std::min_element(main_end.begin(), main_end.end());
        main_f->shift_end(min_end);
        i = 0;
        BOOST_FOREACH (Fragment* o_f, *block) {
            if (o_f != main_f) {
                int delta = main_end[i] - min_end;
                o_f->shift_end(o_end[i] - delta);
                i++;
            }
        }
        const float MIN_ACCEPTED = 0.5;
        if (min_end < batch* MIN_ACCEPTED || min_end == 0) {
            break;
        }
    }
    return result;
}

}

