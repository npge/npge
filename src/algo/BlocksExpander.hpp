/*
 * NPG-explorer, Nucleotide PanGenome explorer
 * Copyright (C) 2014 Boris Nagaev
 *
 * See the LICENSE file for terms of use.
 */

#ifndef NPGE_BLOCKS_EXPANDER_HPP_
#define NPGE_BLOCKS_EXPANDER_HPP_

#include "Processor.hpp"

namespace npge {

/** Add new fragments to blocks.
This adds to the block new fragments, made from neighbor blocks,
if they are \ref PairAligner::aligned() "aligned" with
some fragment from this block.

\warning
   Fragments must be \ref Connector "connected"
   for this to work correctly.
*/
class BlocksExpander : public Processor {
public:
    /** Constructor */
    BlocksExpander();

    /** Expand one block */
    bool expand(Block* block) const;

protected:
    void run_impl() const;
    const char* name_impl() const;
};

}

#endif

