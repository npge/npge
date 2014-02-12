/*
 * bloomrepeats, Find genomic repeats, using Bloom filter based prefiltration
 * Copyright (C) 2012 Boris Nagaev
 *
 * See the LICENSE file for terms of use.
 */

#ifndef BR_DECONSEQ_HPP_
#define BR_DECONSEQ_HPP_

#include "Processor.hpp"

namespace bloomrepeats {

/** Build blocks from blocks, set to sequences of blocks of OtherBlockSet.
Each sequence of OtherBlockSet must have block set.
Those blocks should be aligned.
Blocks of OtherBlockSet are sliced from blocks of sequences and
added to this block set.

\warning Missing sequences are now added to this block.
    They should be added using SequencesFromOther.
*/
class DeConSeq : public Processor {
public:
    /** Constructor */
    DeConSeq(const BlockSetPtr& source = BlockSetPtr());

    /** Produce one block.
    Returned block is not inserted into any block set.
    */
    static Block* deconseq_block(const Block* block);

protected:
    bool run_impl() const;

    const char* name_impl() const;
};

}

#endif

