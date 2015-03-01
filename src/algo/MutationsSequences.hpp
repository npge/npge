/*
 * NPG-explorer, Nucleotide PanGenome explorer
 * Copyright (C) 2012-2015 Boris Nagaev
 *
 * See the LICENSE file for terms of use.
 */

#ifndef NPGE_MUTATIONS_SEQUENCES_HPP_
#define NPGE_MUTATIONS_SEQUENCES_HPP_

#include "BlocksJobs.hpp"
#include "global.hpp"

namespace npge {

class PrintMutations;

/** Create a sequence per genome of variable parts and map blocks */
class MutationsSequences : public BlocksJobs {
public:
    /** Constructor */
    MutationsSequences();

protected:
    ThreadData* before_thread_impl() const;

    void process_block_impl(Block* block, ThreadData* data) const;

    void after_thread_impl(ThreadData* data) const;

    const char* name_impl() const;

private:
    PrintMutations* print_mutations_;
};

}

#endif

