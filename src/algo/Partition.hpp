/*
 * bloomrepeats, Find genomic repeats, using Bloom filter based prefiltration
 * Copyright (C) 2012 Boris Nagaev
 *
 * See the LICENSE file for terms of use.
 */

#ifndef BR_PARTITION_HPP_
#define BR_PARTITION_HPP_

#include "BlocksJobs.hpp"

namespace bloomrepeats {

/** Split fragments of target block set according to fragments of other.

Each new fragment is guaranteed to be a subfragment of a fragment from other
block set.
Ori of new fragments is equal to ori of original fragment from target block set.
*/
class Partition : public BlocksJobs {
public:
    /** Constructor */
    Partition();

    /** Destructor */
    ~Partition();

protected:
    bool change_blocks_impl(std::vector<Block*>& blocks) const;

    bool apply_to_block_impl(Block* block) const;

    const char* name_impl() const;

private:
    struct Impl;

    Impl* impl_; // noncopyable because Processor is noncopyable
};

}

#endif

