/*
 * bloomrepeats, Find genomic repeats, using Bloom filter based prefiltration
 * Copyright (C) 2014 Boris Nagaev
 *
 * See the LICENSE file for terms of use.
 */

#ifndef BR_MULTIPLE_ALIGNER_HPP_
#define BR_MULTIPLE_ALIGNER_HPP_

#include "AbstractAligner.hpp"
#include "global.hpp"

namespace bloomrepeats {

/** Align blocks with internal alignment algorithm */
class MultipleAligner : public AbstractAligner {
public:
    /** Constructor */
    MultipleAligner();

    /** Align multiple sequences */
    static void multiple_aligner(Strings& seqs, PairAligner* pa);

protected:
    const char* name_impl() const;

    void align_seqs_impl(Strings& seqs) const;
};

}

#endif

