/*
 * bloomrepeats, Find genomic repeats, using Bloom filter based prefiltration
 * Copyright (C) 2012 Boris Nagaev
 *
 * See the LICENSE file for terms of use.
 */

#ifndef BR_FRAGMENTS_EXPANDER_HPP_
#define BR_FRAGMENTS_EXPANDER_HPP_

#include "Processor.hpp"
#include "PairAligner.hpp"

namespace bloomrepeats {

/** Expand all blocks (starting from blocks of large number of fragments) */
class FragmentsExpander : public Processor {
public:
    /** Constructor
    \param batch Length of piece, passed to PairAligner at a time.
    \param ori Direction of expansion. 0 means both.
    \param max_overlap Max number of positions, that are allowed to be added
       to the block after first overlap occured.
       -1 means "overlaps of any length are allowed".
       Fragments must be \ref Connector "connected"
       for this to work correctly.

    Steps:
     - One fragment is selected as main.
     - On each iteration, other fragments are aligned to main one.
     - If at least one fragment was aligned on less then 0.5 of batch,
       expansion is stopped.
    */
    FragmentsExpander(int batch = 100, int ori = 0, int max_overlap = 0);

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

    /** Get direction of expansion */
    int ori() const {
        return ori_;
    }

    /** Set direction of expansion */
    void set_ori(int ori) {
        ori_ = ori;
    }

    /** Get max number of positions added after first overlap occured */
    int max_overlap() const {
        return max_overlap_;
    }

    /** Set max number of positions added after first overlap occured */
    void set_max_overlap(int max_overlap) {
        max_overlap_ = max_overlap;
    }

    /** Expand one block */
    bool expand(Block* block) const;

protected:
    /** Apply the action */
    bool run_impl() const;

private:
    mutable PairAligner aligner_;
    int batch_;
    int ori_;
    int max_overlap_;

    bool expand_end(Block* block) const;
};

}

#endif

