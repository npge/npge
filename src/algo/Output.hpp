/*
 * NPG-explorer, Nucleotide PanGenome explorer
 * Copyright (C) 2014 Boris Nagaev
 *
 * See the LICENSE file for terms of use.
 */

#ifndef NPGE_OUTPUT_HPP_
#define NPGE_OUTPUT_HPP_

#include "AbstractOutput.hpp"

namespace npge {

/** Print blocks in fasta format to file or to stdout */
class Output : public AbstractOutput {
public:
    /** Constructor */
    Output(const std::string& prefix = "out-");

protected:
    const char* name_impl() const;

    void print_block(std::ostream& o, Block* block) const;

    void print_header(std::ostream& o) const;
};

}

#endif

