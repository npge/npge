/*
 * bloomrepeats, Find genomic repeats, using Bloom filter based prefiltration
 * Copyright (C) 2014 Boris Nagaev
 *
 * See the LICENSE file for terms of use.
 */

#include <boost/foreach.hpp>

#include "RemoveWeak.hpp"
#include "BlockSet.hpp"
#include "Block.hpp"
#include "global.hpp"

namespace bloomrepeats {

void RemoveWeak::run_impl() const {
    Blocks blocks(block_set()->begin(), block_set()->end());
    BOOST_FOREACH (Block* block, blocks) {
        if (block->weak()) {
            block_set()->erase(block);
        }
    }
}

}

