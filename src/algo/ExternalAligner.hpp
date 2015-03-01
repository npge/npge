/*
 * NPG-explorer, Nucleotide PanGenome explorer
 * Copyright (C) 2012-2015 Boris Nagaev
 *
 * See the LICENSE file for terms of use.
 */

#ifndef NPGE_EXTERNAL_ALIGNER_HPP_
#define NPGE_EXTERNAL_ALIGNER_HPP_

#include "AbstractAligner.hpp"

namespace npge {

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
    std::string aligner_type() const;

    const char* name_impl() const;

    void align_seqs_impl(Strings& seqs) const;
};

/** Mafft aligner */
class MafftAligner : public ExternalAligner {
public:
    /** Constructor */
    MafftAligner();

protected:
    std::string aligner_type() const;

    const char* name_impl() const;
};

/** Muscle aligner */
class MuscleAligner : public ExternalAligner {
public:
    /** Constructor */
    MuscleAligner();

protected:
    std::string aligner_type() const;

    const char* name_impl() const;
};

}

#endif

