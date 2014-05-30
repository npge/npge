/*
 * NPG-explorer, Nucleotide PanGenome explorer
 * Copyright (C) 2014 Boris Nagaev
 *
 * See the LICENSE file for terms of use.
 */

#ifndef BR_META_ALIGNER_HPP_
#define BR_META_ALIGNER_HPP_

#include "AbstractAligner.hpp"

namespace npge {

class ExternalAligner;
class MultipleAligner;
class SimilarAligner;
class DummyAligner;

/** Align blocks with one of possible aligners */
class MetaAligner : public AbstractAligner {
public:
    /** Constructor */
    MetaAligner();

protected:
    const char* name_impl() const;

    void align_seqs_impl(Strings& seqs) const;

private:
    ExternalAligner* external_;
    MultipleAligner* multiple_;
    SimilarAligner* similar_;
    DummyAligner* dummy_;
    mutable AbstractAligner* aligner_;

    bool check_type(std::string& m) const;
};

}

#endif

