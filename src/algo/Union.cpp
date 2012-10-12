/*
 * bloomrepeats, Find genomic repeats, using Bloom filter based prefiltration
 * Copyright (C) 2012 Boris Nagaev
 *
 * See the LICENSE file for terms of use.
 */

#include <boost/foreach.hpp>

#include "Union.hpp"
#include "Block.hpp"
#include "Fragment.hpp"
#include "BlockSet.hpp"

namespace bloomrepeats {

Union::Union(const BlockSetPtr& source):
    OtherBlockSet(source)
{ }

Block* Union::clone_block(Block* source) {
    Block* result = new Block(source->name());
    BOOST_FOREACH (Fragment* f, *source) {
        result->insert(new Fragment(*f));
    }
    return result;
}

BlockSetPtr Union::clone_block_set(BlockSetPtr block_set) {
    BlockSetPtr result = boost::make_shared<BlockSet>();
    Union cloner(block_set);
    cloner.apply(result);
    return result;
}

bool Union::run_impl() const {
    BOOST_FOREACH (Block* block, *other()) {
        block_set()->insert(clone_block(block));
    }
    return !other()->empty();
}

}

