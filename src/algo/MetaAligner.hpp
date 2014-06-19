/*
 * NPG-explorer, Nucleotide PanGenome explorer
 * Copyright (C) 2014 Boris Nagaev
 *
 * See the LICENSE file for terms of use.
 */

#ifndef NPGE_META_ALIGNER_HPP_
#define NPGE_META_ALIGNER_HPP_

#include "AbstractAligner.hpp"

namespace npge {

class ExternalAligner;
class MafftAligner;
class MuscleAligner;
class MultipleAligner;
class SimilarAligner;
class DummyAligner;

/** Align blocks with one of possible aligners */
class MetaAligner : public AbstractAligner {
public:
    /** Constructor */
    MetaAligner();

protected:
    std::string aligner_type() const;

    const char* name_impl() const;

    void align_seqs_impl(Strings& seqs) const;

private:
    ExternalAligner* external_;
    MafftAligner* mafft_;
    MuscleAligner* muscle_;
    MultipleAligner* multiple_;
    SimilarAligner* similar_;
    DummyAligner* dummy_;
    mutable AbstractAligner* aligner_;

    bool check_type(std::string& m) const;
};

}

#endif

