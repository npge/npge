/*
 * NPG-explorer, Nucleotide PanGenome explorer
 * Copyright (C) 2014 Boris Nagaev
 *
 * See the LICENSE file for terms of use.
 */

#include <boost/foreach.hpp>

#include "Clear.hpp"
#include "BlockSet.hpp"

namespace npge {

Clear::Clear() {
    add_opt("clear-blocks", "Remove blocks", true);
    add_opt("clear-seqs", "Remove sequences (warning: this makes "
            "dangling pointers from orphaned blocks, if they remain)",
            false);
    add_opt("clear-bsas", "Remove blockset alignments", true);
    declare_bs("target", "Target blockset");
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

