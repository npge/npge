/*
 * bloomrepeats, Find genomic repeats, using Bloom filter based prefiltration
 * Copyright (C) 2012 Boris Nagaev
 *
 * See the LICENSE file for terms of use.
 */

#include <boost/foreach.hpp>
#include <boost/assert.hpp>

#include "BlockSet.hpp"
#include "Block.hpp"

namespace bloomrepeats {

void BlockSet::insert(BlockPtr block) {
#ifndef NDEBUG
    BOOST_FOREACH (BlockPtr b, *this) {
        BOOST_ASSERT(block != b);
    }
#endif
    blocks_.insert(block);
}

void BlockSet::erase(BlockPtr block) {
    blocks_.erase(block);
}

size_t BlockSet::size() const {
    return blocks_.size();
}

bool BlockSet::empty() const {
    return blocks_.empty();
}

bool BlockSet::has(BlockPtr block) const {
    return blocks_.find(block) != blocks_.end();
}

void BlockSet::clear() {
    blocks_.clear();
}

BlockPtr BlockSet::front() const {
    return empty() ? BlockPtr() : *(begin());
}

BlockSet::iterator BlockSet::begin() {
    return blocks_.begin();
}

BlockSet::const_iterator BlockSet::begin() const {
    return blocks_.begin();
}

BlockSet::iterator BlockSet::end() {
    return blocks_.end();
}

BlockSet::const_iterator BlockSet::end() const {
    return blocks_.end();
}

std::ostream& operator<<(std::ostream& o, const BlockSet& block_set) {
    BOOST_FOREACH (BlockPtr block, block_set) {
        o << *block << std::endl;
    }
    return o;
}

}

