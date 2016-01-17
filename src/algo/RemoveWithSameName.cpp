/*
 * NPG-explorer, Nucleotide PanGenome explorer
 * Copyright (C) 2012-2016 Boris Nagaev
 *
 * See the LICENSE file for terms of use.
 */

#include <boost/foreach.hpp>
#include <boost/cast.hpp>

#include "RemoveWithSameName.hpp"
#include "SortedVector.hpp"
#include "BlockSet.hpp"
#include "Block.hpp"

namespace npge {

RemoveWithSameName::RemoveWithSameName() {
    declare_bs("other", "Reference blocks");
    declare_bs("target", "Modified blockset");
}

typedef SortedVector<std::string> SortedStrings;

void RemoveWithSameName::initialize_work_impl() const {
    names_.clear();
    BOOST_FOREACH (Block* b, *other()) {
        names_.push_back(b->name());
    }
    SortedStrings& nn = static_cast<SortedStrings&>(names_);
    nn.sort_unique();
}

struct RWSNData: public ThreadData {
    Blocks removed_;
};

ThreadData* RemoveWithSameName::before_thread_impl() const {
    return new RWSNData;
}

void RemoveWithSameName::process_block_impl(Block* block,
        ThreadData* d) const {
    RWSNData* data = boost::polymorphic_downcast<RWSNData*>(d);
    Blocks& removed = data->removed_;
    SortedStrings& nn = static_cast<SortedStrings&>(names_);
    if (nn.has_elem(block->name())) {
        removed.push_back(block);
    }
}

void RemoveWithSameName::after_thread_impl(ThreadData* d) const {
    RWSNData* data = boost::polymorphic_downcast<RWSNData*>(d);
    Blocks& removed = data->removed_;
    BlockSet& target = *block_set();
    BOOST_FOREACH (Block* b, removed) {
        target.erase(b);
    }
}

void RemoveWithSameName::finish_work_impl() const {
    names_.clear();
}

const char* RemoveWithSameName::name_impl() const {
    return "Remove from target blocks with names from other";
}

}

