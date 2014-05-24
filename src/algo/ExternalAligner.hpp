/*
 * bloomrepeats, Find genomic repeats, using Bloom filter based prefiltration
 * Copyright (C) 2012 Boris Nagaev
 *
 * See the LICENSE file for terms of use.
 */

#ifndef BR_EXTERNAL_ALIGNER_HPP_
#define BR_EXTERNAL_ALIGNER_HPP_

#include "AbstractAligner.hpp"

namespace bloomrepeats {

/** Align blocks with external alignment tool */
class ExternalAligner : public AbstractAligner {
public:
    /** Constructor */
    ExternalAligner();

    /** Apply external aligner to file */
    void align_file(const std::string& input,
                    const std::string& output) const;

    /** Return list of alignment rows from fasta file */
    void read_alignment(Strings& rows,
                        const std::string& file) const;

protected:
    const char* name_impl() const;

    void align_seqs_impl(Strings& seqs) const;
};

}

#endif

