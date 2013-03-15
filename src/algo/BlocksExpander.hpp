/*
 * bloomrepeats, Find genomic repeats, using Bloom filter based prefiltration
 * Copyright (C) 2012 Boris Nagaev
 *
 * See the LICENSE file for terms of use.
 */

#ifndef BR_BLOCKS_EXPANDER_HPP_
#define BR_BLOCKS_EXPANDER_HPP_

#include "Processor.hpp"
#include "ExpanderBase.hpp"

namespace bloomrepeats {

/** Add new fragments to blocks.
This adds to the block new fragments, made from neighbor blocks,
if they are \ref PairAligner::aligned() "aligned" with
some fragment from this block.

\warning
   Fragments must be \ref Connector "connected"
   for this to work correctly.
*/
class BlocksExpander : public Processor, public ExpanderBase {
public:
    /** Constructor
    \param batch Length of piece, passed to PairAligner at a time.
    */
    BlocksExpander(int batch = 100);

    /** Expand one block */
    bool expand(Block* block) const;

protected:
    void add_options_impl(po::options_description& desc) const;

    void apply_options_impl(const po::variables_map& vm);

    /** Apply the action */
    bool run_impl() const;
};

}

#endif

