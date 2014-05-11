/*
 * bloomrepeats, Find genomic repeats, using Bloom filter based prefiltration
 * Copyright (C) 2014 Boris Nagaev
 *
 * See the LICENSE file for terms of use.
 */

#ifndef BR_MOVE_UNCHANGED_HPP_
#define BR_MOVE_UNCHANGED_HPP_

#include <stdint.h> // for uint32_t

#include "BlocksJobs.hpp"
#include "SortedVector.hpp"

namespace bloomrepeats {

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
    mutable SortedVector<uint32_t> hashes_;
};

}

#endif

