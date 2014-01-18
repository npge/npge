/*
 * bloomrepeats, Find genomic repeats, using Bloom filter based prefiltration
 * Copyright (C) 2012 Boris Nagaev
 *
 * See the LICENSE file for terms of use.
 */

#ifndef BR_SUBSTRACT_HPP_
#define BR_SUBSTRACT_HPP_

#include "BlocksJobs.hpp"

namespace bloomrepeats {

/** Remove from target fragments that have overlaps with other */
class Subtract : public BlocksJobs {
public:
    /** Constructor */
    Subtract();

    /** Destructor */
    ~Subtract();

protected:
    bool change_blocks_impl(std::vector<Block*>& blocks) const;

    bool process_block_impl(Block* block, ThreadData*) const;

private:
    struct Impl;

    Impl* impl_; // noncopyable because Processor is noncopyable
};

}

#endif

