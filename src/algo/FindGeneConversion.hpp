/*
 * bloomrepeats, Find genomic repeats, using Bloom filter based prefiltration
 * Copyright (C) 2014 Boris Nagaev
 *
 * See the LICENSE file for terms of use.
 */

#ifndef BR_FIND_GENE_CONVERSION_HPP_
#define BR_FIND_GENE_CONVERSION_HPP_

#include "BlocksJobs.hpp"
#include "global.hpp"

namespace bloomrepeats {

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

    bool change_blocks_impl(std::vector<Block*>& blocks) const;
    ThreadData* before_thread_impl() const;
    bool process_block_impl(Block* block, ThreadData* data) const;
    bool after_thread_impl(ThreadData* data) const;
};

}

#endif

