/*
 * NPG-explorer, Nucleotide PanGenome explorer
 * Copyright (C) 2012-2015 Boris Nagaev
 *
 * See the LICENSE file for terms of use.
 */

#ifndef NPGE_REMOVE_ALIGNMENT_HPP_
#define NPGE_REMOVE_ALIGNMENT_HPP_

#include "BlocksJobs.hpp"

namespace npge {

/** Remove alignment rows of all fragments */
class RemoveAlignment : public BlocksJobs {
public:
    /** Constructor */
    RemoveAlignment();

protected:
    void process_block_impl(Block* block, ThreadData*) const;

    const char* name_impl() const;
};

}

#endif

