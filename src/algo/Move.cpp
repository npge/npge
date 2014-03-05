/*
 * bloomrepeats, Find genomic repeats, using Bloom filter based prefiltration
 * Copyright (C) 2014 Boris Nagaev
 *
 * See the LICENSE file for terms of use.
 */

#include <vector>
#include <boost/foreach.hpp>

#include "Move.hpp"
#include "BlockSet.hpp"

namespace bloomrepeats {

bool Move::run_impl() const {
    BlockSet& o = *other();
    BlockSet& t = *block_set();
    std::vector<Block*> blocks(o.begin(), o.end());
    BOOST_FOREACH (Block* block, blocks) {
        o.detach(block);
        t.insert(block);
    }
    return true; // TODO
}

}

