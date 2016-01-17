/*
 * NPG-explorer, Nucleotide PanGenome explorer
 * Copyright (C) 2012-2016 Boris Nagaev
 *
 * See the LICENSE file for terms of use.
 */

#include <vector>
#include <boost/foreach.hpp>

#include "Move.hpp"
#include "BlockSet.hpp"

namespace npge {

Move::Move() {
    declare_bs("other", "Source from where blocks are moved");
    declare_bs("target", "Destination to where blocks are moved");
}

void Move::run_impl() const {
    BlockSet& o = *other();
    BlockSet& t = *block_set();
    std::vector<Block*> blocks(o.begin(), o.end());
    BOOST_FOREACH (Block* block, blocks) {
        o.detach(block);
        t.insert(block);
    }
}

const char* Move::name_impl() const {
    return "Move all blocks from other blockset to target blockset";
}

}

