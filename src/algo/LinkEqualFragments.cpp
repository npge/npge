/*
 * NPG-explorer, Nucleotide PanGenome explorer
 * Copyright (C) 2012-2016 Boris Nagaev
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

namespace npge {

struct LinkEqualFragments::Impl {
    typedef std::map<Fragment, Fragment*> F2F;
    F2F f2f_;
};

LinkEqualFragments::LinkEqualFragments() {
    impl_ = new Impl;
    declare_bs("target", "Blockset where fragments, equal to ones "
               "from other blockset, are replaced and blocks are "
               "marked as weak");
    declare_bs("other", "Blockset fragments of which are "
               "searched in target blockset");
}

LinkEqualFragments::~LinkEqualFragments() {
    delete impl_;
    impl_ = 0;
}

void LinkEqualFragments::change_blocks_impl(std::vector<Block*>&) const {
    BOOST_FOREACH (Block* b, *other()) {
        BOOST_FOREACH (Fragment* f, *b) {
            impl_->f2f_[*f] = f;
        }
    }
}

void LinkEqualFragments::process_block_impl(Block* block, ThreadData*) const {
    // block from target
    if (block->empty() || block->weak()) {
        return;
    }
    std::vector<Fragment*> copies;
    BOOST_FOREACH (Fragment* fragment, *block) {
        Impl::F2F::const_iterator it = impl_->f2f_.find(*fragment);
        if (it == impl_->f2f_.end()) {
            // one of fragments can't be replaced
            return;
        }
        copies.push_back(it->second);
    }
    block->clear();
    block->set_weak(true);
    BOOST_FOREACH (Fragment* fragment, copies) {
        block->insert(fragment);
    }
}

const char* LinkEqualFragments::name_impl() const {
    return "Link fragments equal to fragments from other blockset";
}

}

