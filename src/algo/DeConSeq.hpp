/*
 * NPG-explorer, Nucleotide PanGenome explorer
 * Copyright (C) 2012-2015 Boris Nagaev
 *
 * See the LICENSE file for terms of use.
 */

#ifndef NPGE_DECONSEQ_HPP_
#define NPGE_DECONSEQ_HPP_

#include "BlocksJobs.hpp"

namespace npge {

/** Build blocks from blocks, set to sequences of blocks of OtherBlockSet.
Each sequence of OtherBlockSet must have blockset.
Those blocks should be aligned.
Blocks of OtherBlockSet are sliced from blocks of sequences and
added to this blockset.

\warning Missing sequences are now added to this block.
    They should be added using SequencesFromOther.
*/
class DeConSeq : public BlocksJobs {
public:
    /** Constructor */
    DeConSeq(const BlockSetPtr& source = BlockSetPtr());

    /** Produce one block.
    Returned block is not inserted into any blockset.
    */
    static Block* deconseq_block(const Block* block);

protected:
    ThreadData* before_thread_impl() const;

    void process_block_impl(Block* b, ThreadData* d) const;

    void after_thread_impl(ThreadData* d) const;

    const char* name_impl() const;
};

}

#endif

