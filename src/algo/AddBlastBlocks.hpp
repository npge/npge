/*
 * NPG-explorer, Nucleotide PanGenome explorer
 * Copyright (C) 2014 Boris Nagaev
 *
 * See the LICENSE file for terms of use.
 */

#ifndef NPGE_ADD_BLAST_BLOCKS_HPP_
#define NPGE_ADD_BLAST_BLOCKS_HPP_

#include "Pipe.hpp"

namespace npge {

/** Add blast hits as blocks.
Source block set is used as input.
Block set of AddBlastBlocks is used as output.

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

