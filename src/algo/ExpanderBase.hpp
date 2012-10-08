/*
 * bloomrepeats, Find genomic repeats, using Bloom filter based prefiltration
 * Copyright (C) 2012 Boris Nagaev
 *
 * See the LICENSE file for terms of use.
 */

#ifndef BR_EXPANDER_BASE_HPP_
#define BR_EXPANDER_BASE_HPP_

#include "PairAligner.hpp"

namespace bloomrepeats {

/** Base class wrapping PairAligner and batch size */
class ExpanderBase {
public:
    /** Constructor
    \param batch Length of piece, passed to PairAligner at a time.
    */
    ExpanderBase(int batch = 100);

    /** Get pair aligner */
    PairAligner& aligner() const {
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

private:
    mutable PairAligner aligner_;
    int batch_;
};

}

#endif

