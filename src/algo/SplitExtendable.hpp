/*
 * NPG-explorer, Nucleotide PanGenome explorer
 * Copyright (C) 2014 Boris Nagaev
 *
 * See the LICENSE file for terms of use.
 */

#ifndef NPGE_SPLIT_EXTENDABLE_HPP_
#define NPGE_SPLIT_EXTENDABLE_HPP_

#include "BlocksJobs.hpp"

namespace npge {

class FragmentsExtender;
class Filter;

/** Find extendable subblocks.
Find extendable subblocks of original blocks
('other' blockset), add results of extending to target
blockset.
Currently finds blocks of size 2.
*/
class SplitExtendable : public BlocksJobs {
public:
    /** Constructor */
    SplitExtendable();

protected:
    void initialize_work_impl() const;
    ThreadData* before_thread_impl() const;
    void process_block_impl(Block*, ThreadData*) const;
    void after_thread_impl(ThreadData* data) const;
    const char* name_impl() const;

private:
    FragmentsExtender* extender_;
    Filter* filter_;
};

}

#endif

