/*
 * bloomrepeats, Find genomic repeats, using Bloom filter based prefiltration
 * Copyright (C) 2012 Boris Nagaev
 *
 * See the LICENSE file for terms of use.
 */

#ifndef BR_ADD_BLAST_BLOCKS_HPP_
#define BR_ADD_BLAST_BLOCKS_HPP_

#include "Pipe.hpp"

namespace bloomrepeats {

/** Add blast hits as blocks.
Source block set is used as input.
Block set of AddBlastBlocks is used as output.

Blocks must be aligned before this processor.
*/
class AddBlastBlocks : public Pipe {
public:
    /** Constructor */
    AddBlastBlocks(BlockSetPtr source = BlockSetPtr());
};

}

#endif

