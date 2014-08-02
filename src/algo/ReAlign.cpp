/*
 * NPG-explorer, Nucleotide PanGenome explorer
 * Copyright (C) 2014 Boris Nagaev
 *
 * See the LICENSE file for terms of use.
 */

#include <memory>
#include <boost/foreach.hpp>

#include "ReAlign.hpp"
#include "MetaAligner.hpp"
#include "Filter.hpp"
#include "Block.hpp"
#include "BlockSet.hpp"
#include "block_stat.hpp"
#include "block_hash.hpp"
#include "throw_assert.hpp"
#include "cast.hpp"

namespace npge {

ReAlign::ReAlign() {
    aligner_ = new MetaAligner;
    aligner_->set_parent(this);
    filter_ = new Filter;
    filter_->set_parent(this);
    declare_bs("target", "Target blockset");
}

struct RAData : public ThreadData {
    Blocks removed_;
    Blocks inserted_;
};

ThreadData* ReAlign::before_thread_impl() const {
    return new RAData;
}

void ReAlign::process_block_impl(Block* b,
                                 ThreadData* d) const {
    if (has_alignment(b)) {
        std::auto_ptr<Block> copy((b->clone()));
        copy->remove_alignment();
        aligner_->align_block(copy.get());
        if (filter_->is_good_block(copy.get()) &&
                copy->identity() > b->identity()) {
            RAData* data = D_CAST<RAData*>(d);
            data->removed_.push_back(b);
            data->inserted_.push_back(copy.release());
        }
    } else {
        aligner_->align_block(b);
    }
}

void ReAlign::after_thread_impl(ThreadData* d) const {
    RAData* data = D_CAST<RAData*>(d);
    BlockSet& t = *block_set();
    BOOST_FOREACH (Block* b, data->removed_) {
        t.erase(b);
    }
    BOOST_FOREACH (Block* b, data->inserted_) {
        t.insert(b);
    }
}

const char* ReAlign::name_impl() const {
    return "Realign blocks";
}

}

