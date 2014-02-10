/*
 * bloomrepeats, Find genomic repeats, using Bloom filter based prefiltration
 * Copyright (C) 2014 Boris Nagaev
 *
 * See the LICENSE file for terms of use.
 */

#ifndef BR_PRINT_MUTATIONS_HPP_
#define BR_PRINT_MUTATIONS_HPP_

#include "AbstractOutput.hpp"
#include "global.hpp"

namespace bloomrepeats {

class PrintMutations : public AbstractOutput {
public:
    /** Constructor */
    PrintMutations();

    /** Print table block - fr - start_pos - stop_pos - change */
    void print_block(std::ostream& o, Block* block) const;

    void print_header(std::ostream& o) const;
};

}

#endif

