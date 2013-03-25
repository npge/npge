/*
 * bloomrepeats, Find genomic repeats, using Bloom filter based prefiltration
 * Copyright (C) 2012 Boris Nagaev
 *
 * See the LICENSE file for terms of use.
 */

#ifndef BR_MOVE_GAPS_HPP_
#define BR_MOVE_GAPS_HPP_

#include "BlocksJobs.hpp"
#include "RowStorage.hpp"

namespace bloomrepeats {

/** Move terminal letters inside.
Exmaple:
Before: "aaaaa-----a". After: "aaaaaa-----".
Length of tail: 1. Length of gap: 5.
*/
class MoveGaps : public BlocksJobs, public RowStorage {
public:
    /** Constructor */
    MoveGaps(int max_tail = 3, float max_tail_to_gap = 1.0);

    /** Return max length of tail */
    int max_tail() const {
        return max_tail_;
    }

    /** Set max length of tail */
    void set_max_tail(int max_tail) {
        max_tail_ = max_tail;
    }

    /** Return max tail length to gap length ratio */
    float max_tail_to_gap() const {
        return max_tail_to_gap_;
    }

    /** Set max tail length to gap length ratio */
    void set_max_tail_to_gap(float max_tail_to_gap) {
        max_tail_to_gap_ = max_tail_to_gap;
    }

protected:
    void add_options_impl(po::options_description& desc) const;

    void apply_options_impl(const po::variables_map& vm);

    bool apply_to_block_impl(Block* block) const;

    const char* name_impl() const;

private:
    int max_tail_;
    float max_tail_to_gap_;
};

}

#endif

