/*
 * bloomrepeats, Find genomic repeats, using Bloom filter based prefiltration
 * Copyright (C) 2012 Boris Nagaev
 *
 * See the LICENSE file for terms of use.
 */

#include "SizeLimits.hpp"

namespace bloomrepeats {

SizeLimits::SizeLimits(int min_fragment_length, int min_block_size):
    min_fragment_length_(min_fragment_length),
    max_fragment_length_(-1),
    min_block_size_(min_block_size),
    max_block_size_(-1),
    min_spreading_(0.0),
    max_spreading_(0.2),
    min_identity_(0.9),
    max_identity_(1.0),
    min_gaps_(0.0),
    max_gaps_(0.2)
{ }

void SizeLimits::allow_everything() {
    set_min_fragment_length(0);
    set_max_fragment_length(-1);
    set_min_block_size(0);
    set_max_block_size(-1);
    set_min_spreading(0.0);
    set_max_spreading(999999.9);
    set_min_identity(0.0);
    set_max_identity(1.0);
    set_min_gaps(0.0);
    set_max_gaps(1.0);
}

void SizeLimits::add_options_impl(po::options_description& desc) const {
    add_unique_options(desc)
    ("min-fragment", po::value<int>()->default_value(min_fragment_length()),
     "Minimal length of fragments in result")
    ("max-fragment", po::value<int>()->default_value(max_fragment_length()),
     "Maximal length of fragments in result (-1 = everything)")
    ("min-block", po::value<int>()->default_value(min_block_size()),
     "Minimal size of blocks in result")
    ("max-block", po::value<int>()->default_value(max_block_size()),
     "Maximal size of blocks in result (-1 = everything)")
    ("min-spreading", po::value<float>()->default_value(min_spreading()),
     "Minimal fragment length spreading ((max - min) / avg)")
    ("max-spreading", po::value<float>()->default_value(max_spreading()),
     "Maximal fragment length spreading ((max - min) / avg)")
    ("min-identity", po::value<float>()->default_value(min_identity()),
     "Minimal block identity (only if alignment is known, "
     "columns without gaps as 1, columns with gaps as 0.5)")
    ("max-identity", po::value<float>()->default_value(max_identity()),
     "Maximal block identity (only if alignment is known, "
     "columns without gaps as 1, columns with gaps as 0.5)")
    ("min-gaps", po::value<float>()->default_value(min_gaps()),
     "Min gap columns percentage (only if alignment is known)")
    ("max-gaps", po::value<float>()->default_value(max_gaps()),
     "Max gap columns percentage (only if alignment is known)")
   ;
}

void SizeLimits::apply_options_impl(const po::variables_map& vm) {
    if (vm.count("min-fragment")) {
        set_min_fragment_length(vm["min-fragment"].as<int>());
    }
    if (vm.count("max-fragment")) {
        set_max_fragment_length(vm["max-fragment"].as<int>());
    }
    if (vm.count("min-block")) {
        set_min_block_size(vm["min-block"].as<int>());
    }
    if (vm.count("max-block")) {
        set_max_block_size(vm["max-block"].as<int>());
    }
    if (vm.count("min-spreading")) {
        set_min_spreading(vm["min-spreading"].as<float>());
    }
    if (vm.count("max-spreading")) {
        set_max_spreading(vm["max-spreading"].as<float>());
    }
    if (vm.count("min-identity")) {
        set_min_identity(vm["min-identity"].as<float>());
    }
    if (vm.count("max-identity")) {
        set_max_identity(vm["max-identity"].as<float>());
    }
    if (vm.count("min-gaps")) {
        set_min_gaps(vm["min-gaps"].as<float>());
    }
    if (vm.count("max-gaps")) {
        set_max_gaps(vm["max-gaps"].as<float>());
    }
}

}

