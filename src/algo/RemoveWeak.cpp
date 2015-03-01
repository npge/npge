/*
 * NPG-explorer, Nucleotide PanGenome explorer
 * Copyright (C) 2012-2015 Boris Nagaev
 *
 * See the LICENSE file for terms of use.
 */

#include <boost/foreach.hpp>

#include "RemoveWeak.hpp"
#include "BlockSet.hpp"
#include "Block.hpp"
#include "global.hpp"

namespace npge {

RemoveWeak::RemoveWeak() {
    declare_bs("target", "Target blockset");
}

void RemoveWeak::run_impl() const {
    Blocks blocks(block_set()->begin(), block_set()->end());
    BOOST_FOREACH (Block* block, blocks) {
        if (block->weak()) {
            block_set()->erase(block);
        }
    }
}

const char* RemoveWeak::name_impl() const {
    return "Remove weak blocks";
}

}

