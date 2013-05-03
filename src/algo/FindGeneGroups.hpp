/*
 * bloomrepeats, Find genomic repeats, using Bloom filter based prefiltration
 * Copyright (C) 2012 Boris Nagaev
 *
 * See the LICENSE file for terms of use.
 */

#ifndef BR_FIND_GENE_GROUPS_HPP_
#define BR_FIND_GENE_GROUPS_HPP_

#include "BlocksJobs.hpp"

namespace bloomrepeats {

/** Find groups of gene parts according to pangenome.
Block sets:
 - 'target' - block set where gene groups are written to.
   Gene groups are \ref Block::weak() "weak" blocks,
   composed from fragments of blocks from 'genes' block set.
   Names of new blocks are derived from names of blocks of pangenome:
   'block_1', 'block_2', etc.
 - 'genes' (const) - blocks of this block set represent genes.
   Fragments inside block are parts of gene, each of them belongs
   to one fragment from 'pangenome' block set.
 - 'pangenome' (const) - blocks represent similar parts of genomes.
*/
class FindGeneGroups : public BlocksJobs {
public:
    /** Constructor */
    FindGeneGroups();

    /** Destructor */
    ~FindGeneGroups();

protected:
    void change_blocks_impl(std::vector<Block*>& blocks) const;

    bool initialize_thread_impl() const;

    bool apply_to_block_impl(Block* block) const;

    bool finish_thread_impl() const;

    bool finish_work_impl() const;

    const char* name_impl() const;

private:
    struct Impl;

    Impl* impl_; // noncopyable because Processor is noncopyable
};

}

#endif

