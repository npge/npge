/*
 * bloomrepeats, Find genomic repeats, using Bloom filter based prefiltration
 * Copyright (C) 2012 Boris Nagaev
 *
 * See the LICENSE file for terms of use.
 */

#include <vector>
#include <boost/foreach.hpp>

#include "Subtract.hpp"
#include "Connector.hpp"
#include "Rest.hpp"
#include "BlockSet.hpp"
#include "Block.hpp"
#include "Fragment.hpp"
#include "FragmentCollection.hpp"
#include "throw_assert.hpp"
#include "global.hpp"

namespace bloomrepeats {

struct Subtract::Impl {
    typedef std::vector<Fragment> VFragments;
    typedef FragmentCollection<Fragment, VFragments> FC;
    FC fc_;
};

Subtract::Subtract() {
    impl_ = new Impl;
}

Subtract::~Subtract() {
    delete impl_;
    impl_ = 0;
}

void Subtract::change_blocks_impl(std::vector<Block*>& /* blocks */) const {
    Connector c;
    c.apply(other());
    Rest r(other());
    BlockSetPtr rest_of_other = new_bs();
    r.apply(rest_of_other);
    c.apply(rest_of_other);
    Rest r1(rest_of_other);
    BlockSetPtr rest_of_rest_of_other = new_bs();
    r1.apply(rest_of_rest_of_other);
    impl_->fc_.clear();
    impl_->fc_.add_bs(*rest_of_rest_of_other);
    impl_->fc_.prepare();
}

void Subtract::process_block_impl(Block* block, ThreadData*) const {
    std::vector<Fragment*> block_fragments(block->begin(), block->end());
    BOOST_FOREACH (Fragment* fragment, block_fragments) {
        if (impl_->fc_.has_overlap(fragment)) {
            delete fragment;
        }
    }
}

}

