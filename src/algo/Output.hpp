/*
 * bloomrepeats, Find genomic repeats, using Bloom filter based prefiltration
 * Copyright (C) 2012 Boris Nagaev
 *
 * See the LICENSE file for terms of use.
 */

#ifndef BR_OUTPUT_HPP_
#define BR_OUTPUT_HPP_

#include "AbstractOutput.hpp"

namespace bloomrepeats {

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

