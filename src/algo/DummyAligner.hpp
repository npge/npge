/*
 * NPG-explorer, Nucleotide PanGenome explorer
 * Copyright (C) 2014 Boris Nagaev
 *
 * See the LICENSE file for terms of use.
 */

#ifndef NPGE_DUMMY_ALIGNER_HPP_
#define NPGE_DUMMY_ALIGNER_HPP_

#include "AbstractAligner.hpp"

namespace npge {

/** Align by adding gaps to sequences */
class DummyAligner : public AbstractAligner {
public:
    /** Constructor */
    DummyAligner();

protected:
    std::string aligner_type() const;

    const char* name_impl() const;

    void align_seqs_impl(Strings& seqs) const;
};

}

#endif

