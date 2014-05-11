/*
 * bloomrepeats, Find genomic repeats, using Bloom filter based prefiltration
 * Copyright (C) 2014 Boris Nagaev
 *
 * See the LICENSE file for terms of use.
 */

#ifndef BR_REMOVE_WITH_SAME_NAME_HPP_
#define BR_REMOVE_WITH_SAME_NAME_HPP_

#include "BlocksJobs.hpp"
#include "global.hpp"

namespace bloomrepeats {

/** Remove from target blocks with names from other */
class RemoveWithSameName : public BlocksJobs {
public:
    /** Constructor */
    RemoveWithSameName();

protected:
    void initialize_work_impl() const;
    ThreadData* before_thread_impl() const;
    void process_block_impl(Block* block, ThreadData*) const;
    void after_thread_impl(ThreadData* data) const;
    void finish_work_impl() const;
    const char* name_impl() const;

private:
    mutable Strings names_; // sorted
};

}

#endif

