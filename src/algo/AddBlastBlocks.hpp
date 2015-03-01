/*
 * NPG-explorer, Nucleotide PanGenome explorer
 * Copyright (C) 2012-2015 Boris Nagaev
 *
 * See the LICENSE file for terms of use.
 */

#ifndef NPGE_ADD_BLAST_BLOCKS_HPP_
#define NPGE_ADD_BLAST_BLOCKS_HPP_

#include "Pipe.hpp"

namespace npge {

/** Add blast hits as blocks.
Source blockset is used as input.
Blockset of AddBlastBlocks is used as output.

Blocks must be aligned before this processor.
*/
class AddBlastBlocks : public Pipe {
public:
    /** Constructor */
    AddBlastBlocks();

protected:
    const char* name_impl() const;
};

}

#endif

