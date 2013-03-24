/*
 * bloomrepeats, Find genomic repeats, using Bloom filter based prefiltration
 * Copyright (C) 2012 Boris Nagaev
 *
 * See the LICENSE file for terms of use.
 */

#include <boost/foreach.hpp>

#include "FragmentsExpander.hpp"
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

FragmentsExpander::FragmentsExpander(int batch, int ori, int max_overlap):
    ExpanderBase(batch),
    ori_(ori), max_overlap_(max_overlap)
{ }

void FragmentsExpander::add_options_impl(po::options_description& desc) const {
    ExpanderBase::add_options_impl(desc);
    add_unique_options(desc)
    ("max-overlap", po::value<int>()->default_value(max_overlap()),
     "max number of positions added after first overlap occur");
}

void FragmentsExpander::apply_options_impl(const po::variables_map& vm) {
    ExpanderBase::apply_options_impl(vm);
    if (vm.count("max-overlap")) {
        if (vm["max-overlap"].as<int>() < 0) {
            throw Exception("'max-overlap' must be >= 0");
        }
        set_max_overlap(vm["max-overlap"].as<int>());
    }
}

bool FragmentsExpander::run_impl() const {
    std::vector<Block*> bs(block_set()->begin(), block_set()->end());
    std::sort(bs.begin(), bs.end(), block_greater_2);
    bool result = false;
    BOOST_FOREACH (Block* block, bs) {
        result |= expand(block);
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
    bool result = false;
    if (ori() == 1 || ori() == 0) {
        result |= expand_end(block);
    }
    if (ori() == -1 || ori() == 0) {
        block->inverse();
        result |= expand_end(block);
        block->inverse();
    }
    return result;
}

bool FragmentsExpander::expand_end(Block* block) const {
    if (block->size() <= 1) {
        return false;
    }
    std::vector<int> main_end(block->size() - 1), o_end(block->size() - 1);
    Fragment* main_f = block->front();
    bool result = false;
    while (true) {
        int max_shift = block->max_shift_end(max_overlap());
        if (max_shift <= 0) {
            break;
        }
        result = true;
        int shift = std::min(batch(), max_shift);
        std::string main_str = main_f->substr(-1, main_f->length() - 1 + shift);
        aligner().set_first(main_str.c_str(), main_str.size());
        int i = 0;
        BOOST_FOREACH (Fragment* o_f, *block) {
            if (o_f != main_f) {
                std::string o_str = o_f->substr(-1, o_f->length() - 1 + shift);
                aligner().set_second(o_str.c_str(), o_str.size());
                aligner().align(main_end[i], o_end[i]);
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
        if (min_end < batch() * MIN_ACCEPTED || min_end == 0) {
            break;
        }
    }
    return result;
}

}

