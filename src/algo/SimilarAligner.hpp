/*
 * bloomrepeats, Find genomic repeats, using Bloom filter based prefiltration
 * Copyright (C) 2014 Boris Nagaev
 *
 * See the LICENSE file for terms of use.
 */

#ifndef BR_SIMILAR_ALIGNER_PROCESSOR_HPP_
#define BR_SIMILAR_ALIGNER_PROCESSOR_HPP_

#include "AbstractAligner.hpp"
#include "config.hpp"

namespace bloomrepeats {

/** Align blocks with high similarity with internal aligner */
class SimilarAligner : public AbstractAligner {
public:
    /** Constructor */
    SimilarAligner();

    /** Align very similar multiple sequences.
    Align similar equal parts of sequences.
    Jump through 1-point mutations: mismatches and gaps,
    if following mismatch_check or gap_check columns are equal.
    Otherwise find closest aligned_check equal fragment
    and continue from identical part of length aligned_check.
    */
    static void similar_aligner(Strings& seqs,
                                int mismatch_check = MISMATCH_CHECK,
                                int gap_check = GAP_CHECK,
                                int aligned_check = ALIGNED_CHECK);

protected:
    const char* name_impl() const;

    void align_seqs_impl(Strings& seqs) const;
};

}

#endif

