/*
 * NPG-explorer, Nucleotide PanGenome explorer
 * Copyright (C) 2012 Boris Nagaev
 *
 * See the LICENSE file for terms of use.
 */

#ifndef BR_FIND_GENE_GROUPS_HPP_
#define BR_FIND_GENE_GROUPS_HPP_

#include "BlocksJobs.hpp"

namespace npge {

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
    void initialize_work_impl() const;

    ThreadData* before_thread_impl() const;

    void process_block_impl(Block* block, ThreadData* td) const;

    void after_thread_impl(ThreadData* td) const;

    const char* name_impl() const;

private:
    struct Impl;

    Impl* impl_; // noncopyable because Processor is noncopyable
};

}

#endif

