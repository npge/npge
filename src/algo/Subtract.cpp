/*
 * bloomrepeats, Find genomic repeats, using Bloom filter based prefiltration
 * Copyright (C) 2012 Boris Nagaev
 *
 * See the LICENSE file for terms of use.
 */

#include <map>
#include <boost/foreach.hpp>

#include "Subtract.hpp"
#include "Connector.hpp"
#include "Rest.hpp"
#include "BlockSet.hpp"
#include "Block.hpp"
#include "Fragment.hpp"
#include "SortedVector.hpp"
#include "throw_assert.hpp"

namespace bloomrepeats {

typedef SortedVector<Fragment> Fragments;
typedef std::map<Sequence*, Fragments> Seq2Fragments;

struct Subtract::Impl {
    Seq2Fragments seq2fragments_;
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
    impl_->seq2fragments_.clear();
    BOOST_FOREACH (Block* block, *rest_of_rest_of_other) {
        BOOST_FOREACH (Fragment* fragment, *block) {
            Sequence* seq = fragment->seq();
            impl_->seq2fragments_[seq].push_back(*fragment);
        }
    }
    BOOST_FOREACH (Seq2Fragments::value_type& s_and_fs, impl_->seq2fragments_) {
        Fragments& fragments = s_and_fs.second;
        fragments.sort();
    }
}

bool Subtract::apply_to_block_impl(Block* block) const {
    bool result = false;
    std::vector<Fragment*> block_fragments(block->begin(), block->end());
    BOOST_FOREACH (Fragment* fragment, block_fragments) {
        Sequence* seq = fragment->seq();
        Seq2Fragments::const_iterator it = impl_->seq2fragments_.find(seq);
        if (it != impl_->seq2fragments_.end()) {
            const Fragments& fragments = it->second;
            BOOST_ASSERT(!fragments.empty());
            Fragments::const_iterator i2 = fragments.lower_bound(*fragment);
            if (i2 != fragments.end() && i2->common_positions(*fragment)) {
                result = true;
                delete fragment;
            } else if (i2 != fragments.begin()) {
                i2--;
                if (i2->common_positions(*fragment)) {
                    result = true;
                    delete fragment;
                }
            }
        }
    }
    return result;
}

}

