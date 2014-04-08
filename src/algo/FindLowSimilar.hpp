/*
 * bloomrepeats, Find genomic repeats, using Bloom filter based prefiltration
 * Copyright (C) 2014 Boris Nagaev
 *
 * See the LICENSE file for terms of use.
 */

#ifndef BR_FIND_LOW_SIMILAR_HPP_
#define BR_FIND_LOW_SIMILAR_HPP_

#include "BlocksJobs.hpp"

namespace bloomrepeats {

/** Find regions of low similarity in blocks */
class FindLowSimilar : public BlocksJobs {
public:
    /** Constructor */
    FindLowSimilar();

protected:
    ThreadData* before_thread_impl() const;
    void process_block_impl(Block* block, ThreadData*) const;
    void after_thread_impl(ThreadData* data) const;
    const char* name_impl() const;
};

}

#endif

