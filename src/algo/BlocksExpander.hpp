/*
 * bloomrepeats, Find genomic repeats, using Bloom filter based prefiltration
 * Copyright (C) 2012 Boris Nagaev
 *
 * See the LICENSE file for terms of use.
 */

#ifndef BR_BLOCKS_EXPANDER_HPP_
#define BR_BLOCKS_EXPANDER_HPP_

#include "Processor.hpp"
#include "PairAligner.hpp"

namespace bloomrepeats {

/** Add new fragments to blocks.
This adds to the block new fragments, made from neighbor blocks,
if they are \ref PairAligner::aligned() "aligned" with
some fragment from this block.

\warning
   Fragments must be \ref Connector "connected"
   for this to work correctly.
*/
class BlocksExpander : public Processor {
public:
    /** Constructor
    \param batch Length of piece, passed to PairAligner at a time.
    */
    BlocksExpander(int batch = 100);

    /** Get pair aligner */
    const PairAligner& aligner() const {
        return aligner_;
    }

    /** Access pair aligner */
    PairAligner& aligner() {
        return aligner_;
    }

    /** Set pair aligner */
    void set_aligner(const PairAligner& aligner) const {
        aligner_ = aligner;
    }

    /** Get length of piece, passed to PairAligner at a time */
    int batch() const {
        return batch_;
    }

    /** Set length of piece, passed to PairAligner at a time */
    void set_batch(int batch) {
        batch_ = batch;
    }

    /** Expand one block */
    bool expand(Block* block) const;

protected:
    /** Apply the action */
    bool run_impl() const;

private:
    mutable PairAligner aligner_;
    int batch_;
};

}

#endif

