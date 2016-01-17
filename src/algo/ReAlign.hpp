/*
 * NPG-explorer, Nucleotide PanGenome explorer
 * Copyright (C) 2012-2016 Boris Nagaev
 *
 * See the LICENSE file for terms of use.
 */

#ifndef NPGE_REALIGN_HPP
#define NPGE_REALIGN_HPP

#include "BlocksJobs.hpp"

namespace npge {

class MetaAligner;
class Filter;

/** Realign blocks.

A block is replaced with its realigned copy if
the copy is still good (Filter) or original block was bad
and the copy's identity is higher.
If a block has no alignment, then it is aligned.
*/
class ReAlign : public BlocksJobs {
public:
    /** Constructor */
    ReAlign();

protected:
    ThreadData* before_thread_impl() const;
    void process_block_impl(Block* block, ThreadData*) const;
    void after_thread_impl(ThreadData* data) const;
    const char* name_impl() const;

private:
    MetaAligner* aligner_;
    Filter* filter_;
};

}

#endif

