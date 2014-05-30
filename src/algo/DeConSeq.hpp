/*
 * NPG-explorer, Nucleotide PanGenome explorer
 * Copyright (C) 2012 Boris Nagaev
 *
 * See the LICENSE file for terms of use.
 */

#ifndef BR_DECONSEQ_HPP_
#define BR_DECONSEQ_HPP_

#include "BlocksJobs.hpp"

namespace npge {

/** Build blocks from blocks, set to sequences of blocks of OtherBlockSet.
Each sequence of OtherBlockSet must have block set.
Those blocks should be aligned.
Blocks of OtherBlockSet are sliced from blocks of sequences and
added to this block set.

\warning Missing sequences are now added to this block.
    They should be added using SequencesFromOther.
*/
class DeConSeq : public BlocksJobs {
public:
    /** Constructor */
    DeConSeq(const BlockSetPtr& source = BlockSetPtr());

    /** Produce one block.
    Returned block is not inserted into any block set.
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

