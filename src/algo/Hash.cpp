/*
 * bloomrepeats, Find genomic repeats, using Bloom filter based prefiltration
 * Copyright (C) 2013 Boris Nagaev
 *
 * See the LICENSE file for terms of use.
 */

#include <iostream>

#include "Hash.hpp"
#include "FileWriter.hpp"
#include "block_hash.hpp"

namespace bloomrepeats {

Hash::Hash():
    file_writer_(this, "hash-file", "Output file with blockset hash") {
    declare_bs("target", "Target blockset");
}

void Hash::run_impl() const {
    std::ostream& out = file_writer_.output();
    out << blockset_hash(*block_set(), workers()) << "\n";
}

const char* Hash::name_impl() const {
    return "Print hash of blockset";
}

}

