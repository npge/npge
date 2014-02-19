/*
 * bloomrepeats, Find genomic repeats, using Bloom filter based prefiltration
 * Copyright (C) 2012 Boris Nagaev
 *
 * See the LICENSE file for terms of use.
 */

#ifndef BR_FRAGMENTS_EXPANDER_HPP_
#define BR_FRAGMENTS_EXPANDER_HPP_

#include "BlocksJobs.hpp"
#include "config.hpp"

namespace bloomrepeats {

/** Expand all blocks (starting from blocks of large number of fragments).
\note If workers() > 1, then overlaps can happen because of races.

Alignment rows of changed blocks are removed.
*/
class FragmentsExpander : public BlocksJobs {
public:
    /** Constructor
    \param batch Length of piece, passed to PairAligner at a time.
    \param ori Direction of expansion. 0 means both.
    \param max_overlap Max number of positions, that are allowed to be added
       to the block after first overlap occured.
       -1 means "overlaps of any length are allowed".
       Fragments must be \ref Connector "connected"
       for this to work correctly.

    Steps:
     - One fragment is selected as main.
     - On each iteration, other fragments are aligned to main one.
     - If at least one fragment was aligned on less then 0.5 of batch,
       expansion is stopped.
    */
    FragmentsExpander(int batch = EXPANDER_BATCH,
                      int ori = 0, int max_overlap = 0);

    /** Expand one block */
    bool expand(Block* block) const;

protected:
    bool change_blocks_impl(std::vector<Block*>& blocks) const;

    bool process_block_impl(Block* block, ThreadData*) const;

    const char* name_impl() const;

private:
    bool expand_end(Block* block, PairAligner& aligner_copy) const;
};

}

#endif

