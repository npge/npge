/*
 * bloomrepeats, Find genomic repeats, using Bloom filter based prefiltration
 * Copyright (C) 2012 Boris Nagaev
 *
 * See the LICENSE file for terms of use.
 */

#include <boost/foreach.hpp>
#include <boost/cast.hpp>

#include "Filter.hpp"
#include "Fragment.hpp"
#include "Block.hpp"
#include "BlockSet.hpp"
#include "block_stat.hpp"

namespace bloomrepeats {

Filter::Filter(int min_fragment_length, int min_block_size):
    SizeLimits(min_fragment_length, min_block_size)
{ }

bool Filter::is_good_fragment(const Fragment* fragment) const {
    return fragment->valid() && fragment->length() >= min_fragment_length() &&
           (fragment->length() <= max_fragment_length() ||
            max_fragment_length() == -1);
}

bool Filter::filter_block(Block* block) const {
    std::vector<Fragment*> block_copy(block->begin(), block->end());
    bool result = false;
    BOOST_FOREACH (Fragment* fragment, block_copy) {
        if (!is_good_fragment(fragment)) {
            fragment->disconnect();
            block->erase(fragment);
            result = true;
        }
    }
    return result;
}

bool Filter::is_good_block(const Block* block) const {
    if (block->size() < min_block_size()) {
        return false;
    }
    if (block->size() > max_block_size() && max_block_size() != -1) {
        return false;
    }
    AlignmentStat al_stat;
    make_stat(al_stat, block);
    if (al_stat.spreading() < min_spreading()) {
        return false;
    }
    if (al_stat.spreading() > max_spreading()) {
        return false;
    }
    if (al_stat.alignment_rows() == block->size()) {
        float identity = block_identity(al_stat);
        int gaps = al_stat.ident_gap() + al_stat.noident_gap();
        float gaps_p = float(gaps) / al_stat.total();
        if (identity < min_identity() || identity > max_identity()) {
            return false;
        }
        if (gaps_p < min_gaps() || gaps_p > max_gaps()) {
            return false;
        }
    }
    return true;
}

void Filter::add_options_impl(po::options_description& desc) const {
    SizeLimits::add_options_impl(desc);
}

void Filter::apply_options_impl(const po::variables_map& vm) {
    SizeLimits::apply_options_impl(vm);
}

class FilterData : public ThreadData {
public:
    std::vector<Block*> blocks_to_erase;
};

ThreadData* Filter::before_thread_impl() const {
    return new FilterData;
}

bool Filter::process_block_impl(Block* block, ThreadData* d) const {
    filter_block(block);
    if (!is_good_block(block)) {
        FilterData* data = boost::polymorphic_downcast<FilterData*>(d);
        data->blocks_to_erase.push_back(block);
    }
    return true; // TODO
}

bool Filter::after_thread_impl(ThreadData* d) const {
    FilterData* data = boost::polymorphic_downcast<FilterData*>(d);
    BlockSet& target = *block_set();
    BOOST_FOREACH (Block* block, data->blocks_to_erase) {
        target.erase(block);
    }
    return true; // TODO
}

const char* Filter::name_impl() const {
    return "Filter blocks";
}

}

