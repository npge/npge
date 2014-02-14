/*
 * bloomrepeats, Find genomic repeats, using Bloom filter based prefiltration
 * Copyright (C) 2012 Boris Nagaev
 *
 * See the LICENSE file for terms of use.
 */

#ifndef BR_STEM_HPP_
#define BR_STEM_HPP_

#include "BlocksJobs.hpp"
#include "SortedVector.hpp"

namespace bloomrepeats {

/** Filter out blocks not represented in at least one of genomes.
Genome name is obtained by Sequence::genome().
*/
class Stem : public BlocksJobs {
public:
    /** Return if the block is good.
    Depends on calculate_genomes().
    */
    bool is_good_block(const Block* block) const;

    /** Prepare internal genomes list */
    void calculate_genomes() const;

protected:
    bool initialize_work_impl() const;

    ThreadData* before_thread_impl() const;

    bool process_block_impl(Block* block, ThreadData* data) const;

    bool after_thread_impl(ThreadData* data) const;

    const char* name_impl() const;

private:
    mutable SortedVector<std::string> genomes_;
};

}

#endif

