/*
 * bloomrepeats, Find genomic repeats, using Bloom filter based prefiltration
 * Copyright (C) 2012 Boris Nagaev
 *
 * See the LICENSE file for terms of use.
 */

#include <map>
#include <algorithm>
#include <boost/foreach.hpp>

#include "LinkEqualFragments.hpp"
#include "BlockSet.hpp"
#include "Block.hpp"
#include "Fragment.hpp"

namespace bloomrepeats {

struct LinkEqualFragments::Impl {
    typedef std::map<Fragment, Fragment*> F2F;
    F2F f2f_;
};

LinkEqualFragments::LinkEqualFragments() {
    impl_ = new Impl;
}

LinkEqualFragments::~LinkEqualFragments() {
    delete impl_;
    impl_ = 0;
}

void LinkEqualFragments::change_blocks_impl(std::vector<Block*>& ) const {
    BOOST_FOREACH (Block* b, *other()) {
        BOOST_FOREACH (Fragment* f, *b) {
            impl_->f2f_[*f] = f;
        }
    }
    BOOST_FOREACH (Block* b, *block_set()) {
        BOOST_FOREACH (Fragment* f, *b) {
            f->disconnect();
        }
    }
}

bool LinkEqualFragments::apply_to_block_impl(Block* block) const {
    // block from target
    if (block->empty() || block->weak()) {
        return false;
    }
    std::vector<Fragment*> copies;
    BOOST_FOREACH (Fragment* fragment, *block) {
        Impl::F2F::const_iterator it = impl_->f2f_.find(*fragment);
        if (it == impl_->f2f_.end()) {
            // one of fragments can't be replaced
            return false;
        }
        copies.push_back(it->second);
    }
    block->clear();
    block->set_weak(true);
    BOOST_FOREACH (Fragment* fragment, copies) {
        block->insert(fragment);
    }
    return true;
}

const char* LinkEqualFragments::name_impl() const {
    return "Link fragments equal to fragments from other block set";
}

}

