/*
 * bloomrepeats, Find genomic repeats, using Bloom filter based prefiltration
 * Copyright (C) 2012 Boris Nagaev
 *
 * See the LICENSE file for terms of use.
 */

#include <boost/foreach.hpp>

#include "Block.hpp"
#include "Fragment.hpp"

namespace bloomrepeats {

Block::Block()
{ }

void Block::insert(FragmentPtr fragment) {
    fragments_.insert(fragment);
    fragment->block_ = shared_from_this();
}

void Block::erase(FragmentPtr fragment) {
    fragments_.erase(fragment);
    fragment->block_.reset();
}

size_t Block::size() const {
    return fragments_.size();
}

bool Block::empty() const {
    return fragments_.empty();
}

bool Block::has(FragmentPtr fragment) const {
    return fragments_.find(fragment) != fragments_.end();
}

void Block::clear() {
    BOOST_FOREACH (FragmentPtr fragment, *this) {
        erase(fragment);
    }
}

FragmentPtr Block::front() const {
    return empty() ? FragmentPtr() : *(begin());
}

Block::Impl::iterator Block::begin() {
    return fragments_.begin();
}

Block::Impl::const_iterator Block::begin() const {
    return fragments_.begin();
}

Block::Impl::iterator Block::end() {
    return fragments_.end();
}

Block::Impl::const_iterator Block::end() const {
    return fragments_.end();
}

}

