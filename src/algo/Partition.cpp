/*
 * bloomrepeats, Find genomic repeats, using Bloom filter based prefiltration
 * Copyright (C) 2012 Boris Nagaev
 *
 * See the LICENSE file for terms of use.
 */

#include <boost/foreach.hpp>

#include "Partition.hpp"
#include "FragmentCollection.hpp"
#include "Block.hpp"

namespace bloomrepeats {

struct Partition::Impl {
    typedef std::vector<Fragment> Fragments;
    typedef FragmentCollection<Fragment, Fragments> FC;
    FC fc_;
};

Partition::Partition() {
    impl_ = new Impl;
}

Partition::~Partition() {
    delete impl_;
    impl_ = 0;
}

void Partition::change_blocks_impl(std::vector<Block*>& /* blocks */) const {
    impl_->fc_.clear();
    impl_->fc_.add_bs(*other());
    impl_->fc_.prepare();
}

bool Partition::apply_to_block_impl(Block* block) const {
    std::vector<Fragment> new_fragments;
    BOOST_FOREACH (Fragment* fragment, *block) {
        impl_->fc_.find_overlaps(new_fragments, fragment);
    }
    block->clear();
    BOOST_FOREACH (const Fragment& fragment, new_fragments) {
        block->insert(new Fragment(fragment));
    }
    return true;
}

}

