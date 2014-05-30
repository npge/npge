/*
 * NPG-explorer, Nucleotide PanGenome explorer
 * Copyright (C) 2014 Boris Nagaev
 *
 * See the LICENSE file for terms of use.
 */

#ifndef BR_SIMILAR_ALIGNER_PROCESSOR_HPP_
#define BR_SIMILAR_ALIGNER_PROCESSOR_HPP_

#include "AbstractAligner.hpp"

namespace npge {

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
    void similar_aligner(Strings& seqs) const;

protected:
    const char* name_impl() const;

    void align_seqs_impl(Strings& seqs) const;
};

}

#endif

