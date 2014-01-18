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
    typedef std::vector<Fragment*> Fragments;
    typedef FragmentCollection<Fragment*, Fragments> FC;
    FC fc_;
};

Partition::Partition() {
    impl_ = new Impl;
}

Partition::~Partition() {
    delete impl_;
    impl_ = 0;
}

bool Partition::change_blocks_impl(std::vector<Block*>& /* blocks */) const {
    impl_->fc_.clear();
    impl_->fc_.add_bs(*other());
    impl_->fc_.prepare();
    return false;
}

bool Partition::process_block_impl(Block* block, ThreadData*) const {
    std::vector<Fragment*> new_fragments;
    BOOST_FOREACH (Fragment* fragment, *block) {
        std::vector<Fragment> overlaps;
        impl_->fc_.find_overlaps(overlaps, fragment);
        BOOST_FOREACH (const Fragment& overlap, overlaps) {
            Fragment* new_fragment = new Fragment(overlap);
            new_fragment->set_ori(fragment->ori());
            new_fragments.push_back(new_fragment);
        }
    }
    block->clear();
    BOOST_FOREACH (Fragment* fragment, new_fragments) {
        block->insert(fragment);
    }
    return true;
}

const char* Partition::name_impl() const {
    return "Split fragments of target according to fragments of other";
}

}

