/*
 * NPG-explorer, Nucleotide PanGenome explorer
 * Copyright (C) 2014 Boris Nagaev
 *
 * See the LICENSE file for terms of use.
 */

#ifndef NPGE_SPLIT_REPEATS_HPP_
#define NPGE_SPLIT_REPEATS_HPP_

#include "BlocksJobs.hpp"
#include "global.hpp"

namespace npge {

class PrintTree;

/** Find splittable blocks-repeats and split them.
Block is a repeat if it includes fragments from same
sequence.
A repeat is splittable if there is a clade of its tree
such that all fragments of this clade do not form repeat,
and at least one fragment from this clade is on sequence
with repeated fragment. Of cource, clade must
include at least 2 fragments.
Max possible clades are selected as new blocks.

Some fragments may remain non-included in new blocks.

Block sets: other => target (weak blocks).

New blocks get names like "<orig_name>g<n>".

PrintTree is used to build tree.
*/
class SplitRepeats : public BlocksJobs {
public:
    /** Constructor */
    SplitRepeats();

protected:
    const char* name_impl() const;

private:
    PrintTree* tree_;

    ThreadData* before_thread_impl() const;
    void process_block_impl(Block* block, ThreadData* data) const;
    void after_thread_impl(ThreadData* data) const;
};

}

#endif

