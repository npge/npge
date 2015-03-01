/*
 * NPG-explorer, Nucleotide PanGenome explorer
 * Copyright (C) 2012-2015 Boris Nagaev
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
    SBlocks hits_;
};

ThreadData* OverlapFinder::before_thread_impl() const {
    return new OFData;
}

void OverlapFinder::process_block_impl(Block* b,
                                       ThreadData* d) const {
    Fragments ff;
    BOOST_FOREACH (Fragment* f, *b) {
        s2f_.find_overlap_fragments(ff, f);
        if (!ff.empty()) {
            OFData* data = D_CAST<OFData*>(d);
            SBlocks& hits = data->hits_;
            hits.push_back(b);
        }
    }
}

void OverlapFinder::after_thread_impl(ThreadData* d) const {
    OFData* data = D_CAST<OFData*>(d);
    SBlocks& hits = data->hits_;
    hits_.extend(hits);
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

