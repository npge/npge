/*
 * bloomrepeats, Find genomic repeats, using Bloom filter based prefiltration
 * Copyright (C) 2012 Boris Nagaev
 *
 * See the LICENSE file for terms of use.
 */

#ifndef BR_CUT_GAPS_HPP_
#define BR_CUT_GAPS_HPP_

#include "BlocksJobs.hpp"
#include "RowStorage.hpp"

namespace bloomrepeats {

/** Cut longest terminal gap.

Alignment is preserved.

Exmaple:
Before: "-aaaaa-----a-". After: "aaaaa-----a".

If no there is no gapless column, then the block is cleared.
*/
class CutGaps : public BlocksJobs, public RowStorage {
protected:
    void add_options_impl(po::options_description& desc) const;

    void apply_options_impl(const po::variables_map& vm);

    bool apply_to_block_impl(Block* block) const;

    const char* name_impl() const;
};

}

#endif

