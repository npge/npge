/*
 * bloomrepeats, Find genomic repeats, using Bloom filter based prefiltration
 * Copyright (C) 2012 Boris Nagaev
 *
 * See the LICENSE file for terms of use.
 */

#ifndef BR_FILTER_HPP_
#define BR_FILTER_HPP_

#include "BlocksJobs.hpp"
#include "SizeLimits.hpp"

namespace bloomrepeats {

/** Filter out short and invalid fragments.
Fragments are removed (and disconnected).
If block contains too few fragments, it is removed as well
with all its fragments.

\see SizeLimits
*/
class Filter : public BlocksJobs, public SizeLimits {
public:
    /** Constructor */
    Filter(int min_fragment_length = 100, int min_block_size = 2);

    /** Return if fragment is good according to this filter */
    bool is_good_fragment(const Fragment* fragment) const;

    /** Process the block (utility method).
    Return if the block was changed.
    */
    bool filter_block(Block* block) const;

    /** Return if block is good according to this filter.
    Apply filter_block() before.
    */
    bool is_good_block(const Block* block) const;

protected:
    /** Add options to options description */
    void add_options_impl(po::options_description& desc) const;

    /** Apply options from variables map */
    void apply_options_impl(const po::variables_map& vm);

    ThreadData* before_thread_impl() const;

    bool process_block_impl(Block* block, ThreadData* data) const;

    bool after_thread_impl(ThreadData* data) const;

    const char* name_impl() const;
};

}

#endif

