/*
 * NPG-explorer, Nucleotide PanGenome explorer
 * Copyright (C) 2014 Boris Nagaev
 *
 * See the LICENSE file for terms of use.
 */

#ifndef NPGE_META_ALIGNER_HPP_
#define NPGE_META_ALIGNER_HPP_

#include <vector>

#include "AbstractAligner.hpp"

namespace npge {

/** Align blocks with one of possible aligners */
class MetaAligner : public AbstractAligner {
public:
    /** Constructor */
    MetaAligner();

    /** Add aligner.
    Ownership is transferred to MetaAligner.
    Default aligners are pre-added.
    */
    void add_aligner(AbstractAligner* aligner);

protected:
    std::string aligner_type() const;

    const char* name_impl() const;

    void align_seqs_impl(Strings& seqs) const;

private:
    std::vector<AbstractAligner*> aligners_;
    mutable AbstractAligner* aligner_;
    mutable std::string last_aligners_;

    bool check_type(std::string& m) const;
};

}

#endif

