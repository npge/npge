/*
 * bloomrepeats, Find genomic repeats, using Bloom filter based prefiltration
 * Copyright (C) 2014 Boris Nagaev
 *
 * See the LICENSE file for terms of use.
 */

#ifndef BR_SELF_OVERLAPS_RESOLVER_HPP_
#define BR_SELF_OVERLAPS_RESOLVER_HPP_

#include "BlocksJobs.hpp"

namespace bloomrepeats {

/** Resolve self-overlaps */
class SelfOverlapsResolver : public BlocksJobs {
protected:
    void change_blocks_impl(std::vector<Block*>& blocks) const;
    void process_block_impl(Block* block, ThreadData*) const;
};

}

#endif

