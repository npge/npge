/*
 * bloomrepeats, Find genomic repeats, using Bloom filter based prefiltration
 * Copyright (C) 2013 Boris Nagaev
 *
 * See the LICENSE file for terms of use.
 */

#include <iostream>

#include "Hash.hpp"
#include "block_hash.hpp"

namespace bloomrepeats {

Hash::Hash() {
    declare_bs("target", "Target blockset");
}

void Hash::run_impl() const {
    std::cerr << blockset_hash(*block_set(), workers()) << "\n";
}

const char* Hash::name_impl() const {
    return "Print hash of blockset";
}

}

