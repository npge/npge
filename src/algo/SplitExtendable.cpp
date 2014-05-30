/*
 * NPG-explorer, Nucleotide PanGenome explorer
 * Copyright (C) 2014 Boris Nagaev
 *
 * See the LICENSE file for terms of use.
 */

#include <set>
#include <boost/cast.hpp>
#include <boost/scoped_ptr.hpp>
#include <boost/foreach.hpp>

#include "SplitExtendable.hpp"
#include "FragmentsExtender.hpp"
#include "Filter.hpp"
#include "Fragment.hpp"
#include "Block.hpp"
#include "BlockSet.hpp"
#include "throw_assert.hpp"

namespace npge {

SplitExtendable::SplitExtendable() {
    filter_ = new Filter;
    filter_->set_parent(this);
    extender_ = new FragmentsExtender;
    extender_->set_parent(this);
    set_block_set_name("other");
    declare_bs("other", "source blocks");
    declare_bs("target", "output, extended good blocks");
}

void SplitExtendable::initialize_work_impl() const {
    extender_->set_workers(1);
    filter_->set_workers(1);
}

struct SEData: public ThreadData {
    Blocks extended_;
};

ThreadData* SplitExtendable::before_thread_impl() const {
    return new SEData;
}

void try_extend(Fragment* a, Fragment* b,
                Blocks& good_blocks,
                FragmentsExtender* extender,
                Filter* filter) {
    boost::scoped_ptr<Block> block((new Block));
    block->insert(a->clone());
    block->insert(b->clone());
    extender->extend(block.get());
    filter->find_good_subblocks(block.get(), good_blocks);
}

typedef std::set<Fragment*> FragmentSet;

void find_extendable(FragmentSet& ff, Blocks& e,
                     FragmentsExtender* extender,
                     Filter* filter) {
    ASSERT_GTE(ff.size(), 2);
    Fragment* a = *ff.begin();
    ff.erase(a);
    BOOST_FOREACH (Fragment* f, ff) {
        Blocks gb;
        try_extend(a, f, gb, extender, filter);
        if (!gb.empty()) {
            e.insert(e.end(), gb.begin(), gb.end());
            ff.erase(f);
            return;
        }
    }
}

void SplitExtendable::process_block_impl(Block* block,
        ThreadData* d) const {
    FragmentSet ff((block->begin()), block->end());
    SEData* data = boost::polymorphic_downcast<SEData*>(d);
    Blocks& extended = data->extended_;
    while (ff.size() >= 2) {
        find_extendable(ff, extended, extender_, filter_);
    }
}

void SplitExtendable::after_thread_impl(ThreadData* d) const {
    SEData* data = boost::polymorphic_downcast<SEData*>(d);
    BlockSet& target = *block_set();
    BOOST_FOREACH (Block* b, data->extended_) {
        target.insert(b);
    }
}

const char* SplitExtendable::name_impl() const {
    return "Find extendable subblocks";
}

}

