/*
 * NPG-explorer, Nucleotide PanGenome explorer
 * Copyright (C) 2014 Boris Nagaev
 *
 * See the LICENSE file for terms of use.
 */

#ifndef NPGE_SUBSTRACT_HPP_
#define NPGE_SUBSTRACT_HPP_

#include "BlocksJobs.hpp"

namespace npge {

/** Remove from target fragments that have overlaps with other */
class Subtract : public BlocksJobs {
public:
    /** Constructor */
    Subtract();

    /** Destructor */
    ~Subtract();

protected:
    void change_blocks_impl(std::vector<Block*>& blocks) const;

    ThreadData* before_thread_impl() const;

    void process_block_impl(Block* block, ThreadData*) const;

    void after_thread_impl(ThreadData* data) const;

    const char* name_impl() const;

private:
    struct Impl;

    Impl* impl_; // noncopyable because Processor is noncopyable
};

}

#endif

