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

/** Data attached to the thread */
class ThreadData {
public:
    /** Constructor */
    ThreadData();

    /** Destructor */
    virtual ~ThreadData();
};

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
    before running process_block() on them.

    Pre-action.
    */
    void change_blocks(std::vector<Block*>& blocks) const;

    /** Sort blocks by size, length, name */
    void sort_blocks(std::vector<Block*>& blocks) const;

    /** Do something before the work.
    It is applied after change_blocks().

    Pre-action.
    */
    void initialize_work() const;

    /** Do some job before creation of any thread.
    Return pointer to ThreadData which will be passed to
    methods below.
    This can be used to accumulate some data in thread.
    ThreadData* can be 0.
    Pre-action.
    */
    ThreadData* before_thread() const;

    /** Do some job after creation the thread.
    Pre-action.
    */
    void initialize_thread(ThreadData* data) const;

    /** Apply an action to a block */
    void process_block(Block* block, ThreadData* data) const;

    /** Do some job before finish the thread.
    Post-action.

    Does nothing by default.
    */
    void finish_thread(ThreadData* data) const;

    /** Do some job after all threads finished.
    Post-action.

    Does nothing by default.

    Deletes data.
    */
    void after_thread(ThreadData* data) const;

    /** Do some job after applying the action to all blocks.
    Post-action.
    */
    void finish_work() const;

protected:
    void run_impl() const;

    /** Change list of blocks.
    Does nothing by default.
    */
    virtual void change_blocks_impl(std::vector<Block*>& blocks) const;

    /** Do something before other work
    Does nothing by default.
    */
    virtual void initialize_work_impl() const;

    /** Do some job before creation of any thread (implementation).
    Returns 0.
    */
    virtual ThreadData* before_thread_impl() const;

    /** Do some job after creation the thread (implementation).
    Pre-action.

    Does nothing by default.
    */
    virtual void initialize_thread_impl(ThreadData* data) const;

    /** Apply an action to a block (implementation).
    This implementation does nothing.
    */
    virtual void process_block_impl(Block* block, ThreadData* data) const;

    /** Do some job before finish the thread (implementation).
    Post-action.

    Does nothing by default.
    */
    virtual void finish_thread_impl(ThreadData* data) const;

    /** Do some job after all threads finished.
    Does nothing.
    */
    virtual void after_thread_impl(ThreadData* data) const;

    /** Do some job after applying the action to all blocks (implementation).
    Post-action.

    Does nothing by default.
    */
    virtual void finish_work_impl() const;

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

