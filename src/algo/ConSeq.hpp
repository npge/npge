/*
 * bloomrepeats, Find genomic repeats, using Bloom filter based prefiltration
 * Copyright (C) 2012 Boris Nagaev
 *
 * See the LICENSE file for terms of use.
 */

#ifndef BR_CONSEQ_HPP_
#define BR_CONSEQ_HPP_

#include "BlocksJobs.hpp"

namespace bloomrepeats {

/** Add consensus sequences, produced from blocks of source block set.
Depends on UniqueNames. Blocks must be aligned.
*/
class ConSeq : public BlocksJobs {
public:
    /** Constructor */
    ConSeq(const BlockSetPtr& source = BlockSetPtr());

protected:
    ThreadData* before_thread_impl() const;

    void process_block_impl(Block* b, ThreadData* d) const;

    void after_thread_impl(ThreadData* d) const;

    const char* name_impl() const;
};

}

#endif

