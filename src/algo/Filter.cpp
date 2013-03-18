/*
 * bloomrepeats, Find genomic repeats, using Bloom filter based prefiltration
 * Copyright (C) 2012 Boris Nagaev
 *
 * See the LICENSE file for terms of use.
 */

#include <boost/foreach.hpp>

#include "Filter.hpp"
#include "Fragment.hpp"
#include "Block.hpp"
#include "BlockSet.hpp"
#include "block_stat.hpp"

namespace bloomrepeats {

Filter::Filter(int min_fragment_length, int min_block_size):
    SizeLimits(min_fragment_length, min_block_size)
{ }

bool Filter::filter_block(Block* block) const {
    std::vector<Fragment*> block_copy(block->begin(), block->end());
    bool result = false;
    BOOST_FOREACH (Fragment* fragment, block_copy) {
        if (!fragment->valid() || fragment->length() < min_fragment_length()) {
            fragment->disconnect();
            block->erase(fragment);
            result = true;
        }
    }
    return result;
}

void Filter::add_options_impl(po::options_description& desc) const {
    SizeLimits::add_options_impl(desc);
}

void Filter::apply_options_impl(const po::variables_map& vm) {
    SizeLimits::apply_options_impl(vm);
}

bool Filter::run_impl() const {
    bool result = false;
    std::vector<Block*> copy(block_set()->begin(), block_set()->end());
    BOOST_FOREACH (Block* block, copy) {
        result |= filter_block(block);
        bool remove_block = false;
        if (block->size() < min_block_size()) {
            remove_block = true;
        }
        AlignmentStat al_stat;
        make_stat(al_stat, block);
        if (al_stat.spreading > max_spreading()) {
            remove_block = true;
        }
        if (al_stat.alignment_rows == block->size()) {
            float identity = block_identity(al_stat, /* allow gaps */ false);
            int gaps = al_stat.ident_gap + al_stat.noident_gap;
            float gaps_p = float(gaps) / al_stat.total;
            if (identity < min_identity() || gaps_p > max_gaps()) {
                remove_block = true;
            }
        }
        if (remove_block) {
            block_set()->erase(block);
            result = true;
        }
    }
    return result;
}

const char* Filter::name_impl() const {
    return "Filter blocks";
}

}

