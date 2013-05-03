/*
 * bloomrepeats, Find genomic repeats, using Bloom filter based prefiltration
 * Copyright (C) 2012 Boris Nagaev
 *
 * See the LICENSE file for terms of use.
 */

#ifndef BR_BLOCKS_JOBS_HPP_
#define BR_BLOCKS_JOBS_HPP_

#include <vector>

#include "Processor.hpp"

namespace bloomrepeats {

/** Apply an action to each block independently.
Base class.

Blocks should not interfere.
*/
class BlocksJobs : public Processor {
public:
    /** Constructor */
    BlocksJobs(const std::string& block_set_name = "target");

    /** Change list of blocks.
    This action is applied to vist of blocks
    before running apply_to_block() on them.
    */
    void change_blocks(std::vector<Block*>& blocks) const;

    /** Apply an action to a block.
    Return if the block was changed.
    */
    bool apply_to_block(Block* block) const;

protected:
    bool run_impl() const;

    /** Change list of blocks.
    Does nothing by default.
    */
    virtual void change_blocks_impl(std::vector<Block*>& blocks) const;

    /** Apply an action to a block (implementation).
    Return if the block was changed.
    */
    virtual bool apply_to_block_impl(Block* block) const = 0;

    /** Get block set for iteration */
    const std::string& block_set_name() const {
        return block_set_name_;
    }

    /** Set block set for iteration */
    void set_block_set_name(const std::string& block_set_name) {
        block_set_name_ = block_set_name;
    }

private:
    std::string block_set_name_;
};

}

#endif

