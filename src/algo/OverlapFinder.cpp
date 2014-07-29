/*
 * NPG-explorer, Nucleotide PanGenome explorer
 * Copyright (C) 2014 Boris Nagaev
 *
 * See the LICENSE file for terms of use.
 */

#include <boost/foreach.hpp>

#include "OverlapFinder.hpp"
#include "cast.hpp"

namespace npge {

OverlapFinder::OverlapFinder() {
    set_block_set_name("bank");
    declare_bs("bank", "where to find overlapping blocks");
    declare_bs("pattern", "overlaps are looked for with "
               "these blocks");
    declare_bs("hits", "overlapping blocks from bank are "
               "copied here");
}

void OverlapFinder::initialize_work_impl() const {
    s2f_.add_bs(*get_bs("pattern"));
    s2f_.prepare();
}

struct OFData : public ThreadData {
    SortedVector<Block*> hits_;
};

ThreadData* OverlapFinder::before_thread_impl() const {
    return new OFData;
}

void OverlapFinder::process_block_impl(Block* b,
                                       ThreadData* d) const {
    OFData* data = D_CAST<OFData*>(d);
    Fragments ff;
    BOOST_FOREACH (Fragment* f, *b) {
        s2f_.find_overlap_fragments(ff, f);
    }
    SortedVector<Block*>& hits = data->hits_;
    BOOST_FOREACH (Fragment* f, ff) {
        hits.push_back(f->block());
    }
}

void OverlapFinder::after_thread_impl(ThreadData* d) const {
    OFData* data = D_CAST<OFData*>(d);
    hits_.extend(data->hits_);
}

void OverlapFinder::finish_work_impl() const {
    hits_.sort_unique();
    BlockSet& hits = *get_bs("hits");
    BOOST_FOREACH (Block* b, hits_) {
        hits.insert(b->clone());
    }
    hits_.clear();
    s2f_.clear();
}

const char* OverlapFinder::name_impl() const {
    return "Finds blocks from 'bank', "
           "overlapping with 'pattern', "
           "results are copied to 'hits'";
}

}

