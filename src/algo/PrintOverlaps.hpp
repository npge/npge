/*
 * bloomrepeats, Find genomic repeats, using Bloom filter based prefiltration
 * Copyright (C) 2012 Boris Nagaev
 *
 * See the LICENSE file for terms of use.
 */

#ifndef BR_PRINT_OVERLAPS_HPP_
#define BR_PRINT_OVERLAPS_HPP_

#include "AbstractOutput.hpp"

namespace bloomrepeats {

/** Print ASCII diagram with all fragments overlapping with a block.
Fragments must be \ref Connector "connected"

It is recommended to use this processor if blocks have alignment.
*/
class PrintOverlaps : public AbstractOutput {
public:
    /** Constructor */
    PrintOverlaps();

    /** Print ASCII diagram with all fragments overlapping with a block */
    void print_block(std::ostream& o, Block* block) const;

protected:
    const char* name_impl() const;
};

}

#endif

