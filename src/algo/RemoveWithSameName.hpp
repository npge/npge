/*
 * NPG-explorer, Nucleotide PanGenome explorer
 * Copyright (C) 2012-2015 Boris Nagaev
 *
 * See the LICENSE file for terms of use.
 */

#ifndef NPGE_REMOVE_WITH_SAME_NAME_HPP_
#define NPGE_REMOVE_WITH_SAME_NAME_HPP_

#include "BlocksJobs.hpp"
#include "global.hpp"

namespace npge {

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

