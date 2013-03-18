/*
 * bloomrepeats, Find genomic repeats, using Bloom filter based prefiltration
 * Copyright (C) 2012 Boris Nagaev
 *
 * See the LICENSE file for terms of use.
 */

#include "SizeLimits.hpp"

namespace bloomrepeats {

SizeLimits::SizeLimits(int min_fragment_length, int min_block_size):
    min_fragment_length_(min_fragment_length), min_block_size_(min_block_size),
    max_spreading_(0.2), min_identity_(0.9), max_gaps_(0.2)
{ }

void SizeLimits::add_options_impl(po::options_description& desc) const {
    add_unique_options(desc)
    ("min-fragment", po::value<size_t>()->default_value(min_fragment_length()),
     "Minimal length of fragments in result")
    ("min-block", po::value<size_t>()->default_value(min_block_size()),
     "Minimal size of blocks in result")
    ("max-spreading", po::value<float>()->default_value(max_spreading()),
     "Fragment length spreading ((max - min) / avg)")
    ("min-identity", po::value<float>()->default_value(min_identity()),
     "Minimal block identity (on non-gap columns, only if alignment is known)")
    ("max-gaps", po::value<float>()->default_value(max_gaps()),
     "Max gap columns percentage (only if alignment is known)")
   ;
}

void SizeLimits::apply_options_impl(const po::variables_map& vm) {
    set_min_fragment_length(vm["min-fragment"].as<size_t>());
    set_min_block_size(vm["min-block"].as<size_t>());
    set_max_spreading(vm["max-spreading"].as<float>());
    set_min_identity(vm["min-identity"].as<float>());
    set_max_gaps(vm["max-gaps"].as<float>());
}

}

