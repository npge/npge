/*
 * NPG-explorer, Nucleotide PanGenome explorer
 * Copyright (C) 2013 Boris Nagaev
 *
 * See the LICENSE file for terms of use.
 */

#include <iostream>

#include "Hash.hpp"
#include "FileWriter.hpp"
#include "block_hash.hpp"

namespace npge {

Hash::Hash():
    file_writer_(this, "hash-file", "Output file with blockset hash") {
    declare_bs("target", "Target blockset");
}

void Hash::run_impl() const {
    std::ostream& out = file_writer_.output();
    out << blockset_hash(*block_set(), workers()) << std::endl;
}

const char* Hash::name_impl() const {
    return "Print hash of blockset";
}

}

