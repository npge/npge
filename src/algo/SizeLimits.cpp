/*
 * bloomrepeats, Find genomic repeats, using Bloom filter based prefiltration
 * Copyright (C) 2012 Boris Nagaev
 *
 * See the LICENSE file for terms of use.
 */

#include "SizeLimits.hpp"

namespace bloomrepeats {

SizeLimits::SizeLimits(int min_fragment_length, int min_block_size):
    min_fragment_length_(min_fragment_length), min_block_size_(min_block_size)
{ }

void SizeLimits::add_options_impl(po::options_description& desc) const {
    add_unique_options(desc)
    ("min-fragment", po::value<size_t>()->default_value(min_fragment_length()),
     "Minimal length of fragments in result")
    ("min-block", po::value<size_t>()->default_value(min_block_size()),
     "Minimal size of blocks in result")
   ;
}

void SizeLimits::apply_options_impl(const po::variables_map& vm) {
    set_min_fragment_length(vm["min-fragment"].as<size_t>());
    set_min_block_size(vm["min-block"].as<size_t>());
}

}

