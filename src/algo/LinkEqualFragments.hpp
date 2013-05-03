/*
 * bloomrepeats, Find genomic repeats, using Bloom filter based prefiltration
 * Copyright (C) 2012 Boris Nagaev
 *
 * See the LICENSE file for terms of use.
 */

#ifndef BR_LINK_EQUAL_FRAGMENTS_HPP_
#define BR_LINK_EQUAL_FRAGMENTS_HPP_

#include "BlocksJobs.hpp"

namespace bloomrepeats {

/** Replace fragments equal to fragments from other block set with links.
Non-weak block all fragments of which are equal to some fragment
from other block set, is marked as weak,
all its fragments are removed (and deleted)
and replaced with pointers to equal fragments from other block set.

\warning Fragments of target block set must not be used anywhere in weak blocks.
\note Fragments of target block are disconnected.
*/
class LinkEqualFragments : public BlocksJobs {
public:
    /** Constructor */
    LinkEqualFragments();

    /** Destructor */
    ~LinkEqualFragments();

protected:
    void change_blocks_impl(std::vector<Block*>& blocks) const;

    bool apply_to_block_impl(Block* block) const;

    const char* name_impl() const;

private:
    struct Impl;

    Impl* impl_; // noncopyable because Processor is noncopyable
};

}

#endif

