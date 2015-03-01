/*
 * NPG-explorer, Nucleotide PanGenome explorer
 * Copyright (C) 2012-2015 Boris Nagaev
 *
 * See the LICENSE file for terms of use.
 */

#ifndef NPGE_FILTER_HPP_
#define NPGE_FILTER_HPP_

#include "BlocksJobs.hpp"

namespace npge {

/** Filter out small fragments or blocks */
class LiteFilter : public BlocksJobs {
public:
    /** Constructor */
    LiteFilter();

protected:
    ThreadData* before_thread_impl() const;

    void process_block_impl(Block* block,
                            ThreadData* data) const;

    void after_thread_impl(ThreadData* data) const;

    const char* name_impl() const;
};

/** Filter out short and invalid fragments.
Fragments are removed (and disconnected).
If block contains too few fragments, it is removed as well
with all its fragments.

\see add_size_limits_options
*/
class Filter : public BlocksJobs {
public:
    /** Constructor */
    Filter();

    /** Return if fragment is good according to this filter */
    bool is_good_fragment(const Fragment* fragment) const;

    /** Process the block (utility method).
    Return if the block was changed.
    */
    bool filter_block(Block* block) const;

    /** Return if block is good according to this filter.
    Apply filter_block() before.
    Alignment must be good in the beginning and in the end of block.
    */
    bool is_good_block(const Block* block) const;

    /** Slice out good subblocks.
    Ownership of new blocks is transferred to caller.
    */
    void find_good_subblocks(const Block* block,
                             Blocks& good_subblocks) const;

protected:
    ThreadData* before_thread_impl() const;

    void process_block_impl(Block* block, ThreadData* data) const;

    void after_thread_impl(ThreadData* data) const;

    const char* name_impl() const;
};

}

#endif

