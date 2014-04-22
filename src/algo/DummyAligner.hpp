/*
 * bloomrepeats, Find genomic repeats, using Bloom filter based prefiltration
 * Copyright (C) 2012 Boris Nagaev
 *
 * See the LICENSE file for terms of use.
 */

#ifndef BR_DUMMY_ALIGNER_HPP_
#define BR_DUMMY_ALIGNER_HPP_

#include "AbstractAligner.hpp"

namespace bloomrepeats {

/** Align by adding gaps to sequences */
class DummyAligner : public AbstractAligner {
public:
    /** Constructor */
    DummyAligner();

protected:
    const char* name_impl() const;

    void align_seqs_impl(Strings& seqs) const;
};

}

#endif

