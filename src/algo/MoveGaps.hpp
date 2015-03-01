/*
 * NPG-explorer, Nucleotide PanGenome explorer
 * Copyright (C) 2012-2015 Boris Nagaev
 *
 * See the LICENSE file for terms of use.
 */

#ifndef NPGE_MOVE_GAPS_HPP_
#define NPGE_MOVE_GAPS_HPP_

#include "BlocksJobs.hpp"

namespace npge {

/** Move terminal letters inside.
Exmaple:
Before: "aaaaa-----a". After: "aaaaaa-----".
Length of tail: 1. Length of gap: 5.
*/
class MoveGaps : public BlocksJobs {
public:
    /** Constructor */
    MoveGaps();

    /** Do the job and return if the block was changed */
    bool move_gaps(Block* block) const;

protected:
    void process_block_impl(Block* block, ThreadData*) const;

    const char* name_impl() const;
};

}

#endif

