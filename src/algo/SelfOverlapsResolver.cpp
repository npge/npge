/*
 * bloomrepeats, Find genomic repeats, using Bloom filter based prefiltration
 * Copyright (C) 2014 Boris Nagaev
 *
 * See the LICENSE file for terms of use.
 */

#include <boost/foreach.hpp>

#include "SelfOverlapsResolver.hpp"
#include "hit.hpp"

namespace bloomrepeats {

bool SelfOverlapsResolver::change_blocks_impl(
    std::vector<Block*>& blocks) const {
    BOOST_FOREACH (Block* block, blocks) {
        BOOST_FOREACH (Fragment* f, *block) {
            f->disconnect();
        }
    }
}

bool SelfOverlapsResolver::process_block_impl(Block* block,
        ThreadData*) const {
    fix_self_overlaps(block);
}

}

