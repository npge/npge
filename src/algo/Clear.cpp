/*
 * bloomrepeats, Find genomic repeats, using Bloom filter based prefiltration
 * Copyright (C) 2012 Boris Nagaev
 *
 * See the LICENSE file for terms of use.
 */

#include <boost/foreach.hpp>

#include "Clear.hpp"
#include "BlockSet.hpp"

namespace bloomrepeats {

Clear::Clear() {
    add_opt("clear-blocks", "Remove blocks", true);
    add_opt("clear-seqs", "Remove sequences (warning: this makes "
            "dangling pointers from orphaned blocks, if they remain)",
            false);
    add_opt("clear-bsas", "Remove block set alignments", true);
}

void Clear::run_impl() const {
    if (opt_value("clear-blocks").as<bool>()) {
        block_set()->clear_blocks();
    }
    if (opt_value("clear-seqs").as<bool>()) {
        block_set()->clear_seqs();
    }
    if (opt_value("clear-bsas").as<bool>()) {
        block_set()->clear_bsas();
    }
}

const char* Clear::name_impl() const {
    return "Remove all blocks and/or sequences from blockset";
}

}

