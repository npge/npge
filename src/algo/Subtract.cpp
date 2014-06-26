/*
 * NPG-explorer, Nucleotide PanGenome explorer
 * Copyright (C) 2014 Boris Nagaev
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
#include "to_s.hpp"
#include "global.hpp"

namespace npge {

struct Subtract::Impl {
    typedef std::vector<Fragment> VFragments;
    typedef FragmentCollection<Fragment, VFragments> FC;
    FC fc_;
};

Subtract::Subtract() {
    impl_ = new Impl;
    add_opt("subtract-equal", "Delete only equal "
            "fragments", false);
    declare_bs("other", "Blocks to which overlaps are found");
    declare_bs("target", "Blocks which are removed if overlap "
               "with blocks from other");
}

Subtract::~Subtract() {
    delete impl_;
    impl_ = 0;
}

void Subtract::change_blocks_impl(std::vector<Block*>& /* blocks */) const {
    impl_->fc_.clear();
    impl_->fc_.add_bs(*other());
    impl_->fc_.prepare();
}

typedef std::pair<Block*, Fragment*> BF;

struct SData : public ThreadData {
    std::vector<BF> to_erase_;
};

ThreadData* Subtract::before_thread_impl() const {
    return new SData;
}

static bool positions_equal(const Fragment* f1,
                            const Fragment* f2) {
    return f1->min_pos() == f2->min_pos() &&
           f1->max_pos() == f2->max_pos() &&
           f1->seq() == f2->seq();
}

void Subtract::process_block_impl(Block* block,
                                  ThreadData* td) const {
    bool equal = opt_value("subtract-equal").as<bool>();
    SData* sd = D_CAST<SData*>(td);
    Fragments block_fragments(block->begin(), block->end());
    BOOST_FOREACH (Fragment* fragment, block_fragments) {
        if (equal) {
            Fragments oo;
            impl_->fc_.find_overlap_fragments(oo, fragment);
            bool to_delete = false;
            BOOST_FOREACH (Fragment* o, oo) {
                if (positions_equal(o, fragment)) {
                    to_delete = true;
                    break;
                }
            }
            if (to_delete) {
                sd->to_erase_.push_back(BF(block, fragment));
            }
        } else if (impl_->fc_.has_overlap(fragment)) {
            sd->to_erase_.push_back(BF(block, fragment));
        }
    }
}

void Subtract::after_thread_impl(ThreadData* td) const {
    SData* sd = D_CAST<SData*>(td);
    BOOST_FOREACH (const BF& bf, sd->to_erase_) {
        Block* block = bf.first;
        Fragment* fragment = bf.second;
        block->erase(fragment);
    }
}

const char* Subtract::name_impl() const {
    return "Remove from target fragments that have overlaps "
           "with other";
}

}

