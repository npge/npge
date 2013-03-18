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
        if (block->size() < min_block_size()) {
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

