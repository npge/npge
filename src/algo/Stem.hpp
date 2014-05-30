/*
 * NPG-explorer, Nucleotide PanGenome explorer
 * Copyright (C) 2014 Boris Nagaev
 *
 * See the LICENSE file for terms of use.
 */

#ifndef NPGE_STEM_HPP_
#define NPGE_STEM_HPP_

#include "BlocksJobs.hpp"
#include "SortedVector.hpp"

namespace npge {

/** Filter out blocks not represented in at least one of genomes.
Genome name is obtained by Sequence::genome().
*/
class Stem : public BlocksJobs {
public:
    /** Constructor */
    Stem();

    /** Return if the block is good.
    Depends on calculate_genomes().
    */
    bool is_good_block(const Block* block) const;

    /** Prepare internal genomes list */
    void calculate_genomes() const;

protected:
    void initialize_work_impl() const;

    ThreadData* before_thread_impl() const;

    void process_block_impl(Block* block, ThreadData* data) const;

    void after_thread_impl(ThreadData* data) const;

    const char* name_impl() const;

private:
    mutable SortedVector<std::string> genomes_;
    mutable bool exact_;
};

}

#endif

