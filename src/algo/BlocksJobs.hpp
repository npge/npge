/*
 * bloomrepeats, Find genomic repeats, using Bloom filter based prefiltration
 * Copyright (C) 2012 Boris Nagaev
 *
 * See the LICENSE file for terms of use.
 */

#ifndef BR_BLOCKS_JOBS_HPP_
#define BR_BLOCKS_JOBS_HPP_

#include "Processor.hpp"

namespace bloomrepeats {

/** Apply an action to each block independently.
Base class.

Blocks should not interfere.
*/
class BlocksJobs : public Processor {
public:
    /** Apply an action to a block.
    Return if the block was changed.
    */
    bool apply_to_block(Block* block) const;

protected:
    bool run_impl() const;

    /** Apply an action to a block (implementation).
    Return if the block was changed.
    */
    virtual bool apply_to_block_impl(Block* block) const = 0;
};

}

#endif

