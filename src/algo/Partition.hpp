/*
 * NPG-explorer, Nucleotide PanGenome explorer
 * Copyright (C) 2012 Boris Nagaev
 *
 * See the LICENSE file for terms of use.
 */

#ifndef NPGE_PARTITION_HPP_
#define NPGE_PARTITION_HPP_

#include "BlocksJobs.hpp"

namespace npge {

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
    void change_blocks_impl(std::vector<Block*>& blocks) const;

    void process_block_impl(Block* block, ThreadData*) const;

    const char* name_impl() const;

private:
    struct Impl;

    Impl* impl_; // noncopyable because Processor is noncopyable
};

}

#endif

