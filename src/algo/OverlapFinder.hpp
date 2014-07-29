/*
 * NPG-explorer, Nucleotide PanGenome explorer
 * Copyright (C) 2014 Boris Nagaev
 *
 * See the LICENSE file for terms of use.
 */

#ifndef NPGE_OVERLAPS_FINDER_HPP_
#define NPGE_OVERLAPS_FINDER_HPP_

#include "BlocksJobs.hpp"
#include "FragmentCollection.hpp"
#include "SortedVector.hpp"

namespace npge {

typedef SortedVector<Block*> SBlocks;

/** Finds blocks from "bank", overlapping with "pattern".

Results are copied to "hits".
*/
class OverlapFinder : public BlocksJobs {
public:
    /** Constructor */
    OverlapFinder();

protected:
    void initialize_work_impl() const;

    ThreadData* before_thread_impl() const;

    void process_block_impl(Block* block,
                            ThreadData* data) const;

    void after_thread_impl(ThreadData* data) const;

    void finish_work_impl() const;

    const char* name_impl() const;

private:
    typedef FragmentCollection<Fragment*, Fragments> S2F;

    mutable S2F s2f_;
    mutable SBlocks hits_;
};

}

#endif

