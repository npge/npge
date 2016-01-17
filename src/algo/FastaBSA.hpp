/*
 * NPG-explorer, Nucleotide PanGenome explorer
 * Copyright (C) 2012-2016 Boris Nagaev
 *
 * See the LICENSE file for terms of use.
 */

#ifndef NPGE_FASTA_BLOCK_SET_ALIGNMENT_HPP_
#define NPGE_FASTA_BLOCK_SET_ALIGNMENT_HPP_

#include "global.hpp"
#include "Processor.hpp"
#include "FileWriter.hpp"

namespace npge {

/** Print blockset alignment as fasta */
class FastaBSA : public Processor {
public:
    /** Constructor */
    FastaBSA();

protected:
    void run_impl() const;

    const char* name_impl() const;

private:
    FileWriter file_writer_;
};

}

#endif

