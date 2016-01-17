/*
 * NPG-explorer, Nucleotide PanGenome explorer
 * Copyright (C) 2012-2016 Boris Nagaev
 *
 * See the LICENSE file for terms of use.
 */

#ifndef NPGE_SELF_OVERLAPS_RESOLVER_HPP_
#define NPGE_SELF_OVERLAPS_RESOLVER_HPP_

#include "BlocksJobs.hpp"

namespace npge {

/** Resolve self-overlaps */
class SelfOverlapsResolver : public BlocksJobs {
public:
    /** Constructor */
    SelfOverlapsResolver();

protected:
    void process_block_impl(Block* block, ThreadData*) const;
    const char* name_impl() const;
};

}

#endif

