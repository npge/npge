/*
 * bloomrepeats, Find genomic repeats, using Bloom filter based prefiltration
 * Copyright (C) 2014 Boris Nagaev
 *
 * See the LICENSE file for terms of use.
 */

#ifndef BR_MUTATIONS_SEQUENCES_HPP_
#define BR_MUTATIONS_SEQUENCES_HPP_

#include "BlocksJobs.hpp"
#include "global.hpp"

namespace bloomrepeats {

/** Create a sequence per genome of variable parts and map blocks */
class MutationsSequences : public BlocksJobs {
public:
    /** Constructor */
    MutationsSequences();

protected:
    ThreadData* before_thread_impl() const;

    bool process_block_impl(Block* block, ThreadData* data) const;

    bool after_thread_impl(ThreadData* data) const;

    const char* name_impl() const;

private:
    PrintMutations* print_mutations_;
};

}

#endif

