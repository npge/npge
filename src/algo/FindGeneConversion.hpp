/*
 * NPG-explorer, Nucleotide PanGenome explorer
 * Copyright (C) 2014 Boris Nagaev
 *
 * See the LICENSE file for terms of use.
 */

#ifndef NPGE_FIND_GENE_CONVERSION_HPP_
#define NPGE_FIND_GENE_CONVERSION_HPP_

#include "BlocksJobs.hpp"
#include "global.hpp"

namespace npge {

class FragmentDistance;

/** Find groups of blocks with gene conversion in them.
Groups of fragments with probable gene conversion
are found in "other" blockset.
These fragments are added to "target" blockset as weak blocks.

FragmentDistance is used to calculate distance between fragments.
*/
class FindGeneConversion : public BlocksJobs {
public:
    /** Constructor */
    FindGeneConversion();

private:
    FragmentDistance* distance_;

    void change_blocks_impl(std::vector<Block*>& blocks) const;
    ThreadData* before_thread_impl() const;
    void process_block_impl(Block* block, ThreadData* data) const;
    void after_thread_impl(ThreadData* data) const;
};

}

#endif

