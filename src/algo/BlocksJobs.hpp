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

    Return is some of target blocks was changed.

    Pre-action.
    */
    bool change_blocks(std::vector<Block*>& blocks) const;

    /** Do some job after creation the thread.
    Return is some of target blocks was changed.

    Pre-action.
    */
    bool initialize_thread() const;

    /** Do some job before finish the thread.
    Return is some of target blocks was changed.

    Post-action.

    Does nothing by default and return false.
    */
    bool finish_thread() const;

    /** Do some job after applying the action to all blocks.
    Return is some of target blocks was changed.

    Post-action.
    */
    bool finish_work() const;

protected:
    bool run_impl() const;

    /** Change list of blocks.
    Return is some of target blocks was changed.
    Does nothing by default.
    */
    virtual bool change_blocks_impl(std::vector<Block*>& blocks) const;

    /** Do some job after creation the thread (implementation).
    Return is some of target blocks was changed.

    Pre-action.

    Does nothing by default and return false.
    */
    virtual bool initialize_thread_impl() const;

    /** Apply an action to a block (implementation).
    Return if the block was changed.
    This implementation does nothing.
    */
    bool apply_to_block_impl(Block* block) const;

    /** Do some job before finish the thread (implementation).
    Return is some of target blocks was changed.

    Post-action.

    Does nothing by default and return false.
    */
    virtual bool finish_thread_impl() const;

    /** Do some job after applying the action to all blocks (implementation).
    Return is some of target blocks was changed.

    Post-action.

    Does nothing by default and return false.
    */
    virtual bool finish_work_impl() const;

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

