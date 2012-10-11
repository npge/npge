/*
 * bloomrepeats, Find genomic repeats, using Bloom filter based prefiltration
 * Copyright (C) 2012 Boris Nagaev
 *
 * See the LICENSE file for terms of use.
 */

#include <boost/foreach.hpp>

#include "Union.hpp"
#include "Block.hpp"
#include "BlockSet.hpp"

namespace bloomrepeats {

Union::Union(const BlockSetPtr& source):
    source_(source)
{ }

bool Union::run_impl() const {
    BOOST_FOREACH (Block* block, *source()) {
        block_set()->insert(block->clone());
    }
    return !source()->empty();
}

}

