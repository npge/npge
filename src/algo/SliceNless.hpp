/*
 * NPG-explorer, Nucleotide PanGenome explorer
 * Copyright (C) 2012-2015 Boris Nagaev
 *
 * See the LICENSE file for terms of use.
 */

#ifndef NPGE_SLICE_N_LESS_HPP_
#define NPGE_SLICE_N_LESS_HPP_

#include "BlocksJobs.hpp"

namespace npge {

/** Find blocks with long N's and slice them to good parts */
class SliceNless : public BlocksJobs {
public:
    /** Constructor */
    SliceNless();

protected:
    ThreadData* before_thread_impl() const;

    void process_block_impl(Block* block,
                            ThreadData* data) const;

    void after_thread_impl(ThreadData* data) const;

    const char* name_impl() const;
};

}

#endif

