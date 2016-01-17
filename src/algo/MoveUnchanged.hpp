/*
 * NPG-explorer, Nucleotide PanGenome explorer
 * Copyright (C) 2012-2016 Boris Nagaev
 *
 * See the LICENSE file for terms of use.
 */

#ifndef NPGE_MOVE_UNCHANGED_HPP_
#define NPGE_MOVE_UNCHANGED_HPP_

#include "BlocksJobs.hpp"
#include "SortedVector.hpp"
#include "global.hpp"

namespace npge {

/** Move unchanged blocks from other to target.
Memorize hashes of all blocks. Running nex time, moves
blocks with old hashes from other to target.
*/
class MoveUnchanged : public BlocksJobs {
public:
    /** Constructor */
    MoveUnchanged();

    /** Clear internal hash set.
    \todo How to clear internal hash set from script?
    */
    void clear();

protected:
    ThreadData* before_thread_impl() const;

    void process_block_impl(Block* block, ThreadData* data) const;

    void after_thread_impl(ThreadData* data) const;

    void finish_work_impl() const;

    const char* name_impl() const;

private:
    mutable SortedVector<hash_t> hashes_;
};

}

#endif

