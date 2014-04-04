/*
 * bloomrepeats, Find genomic repeats, using Bloom filter based prefiltration
 * Copyright (C) 2012 Boris Nagaev
 *
 * See the LICENSE file for terms of use.
 */

#ifndef BR_ORI_BY_MAJORITY_HPP_
#define BR_ORI_BY_MAJORITY_HPP_

#include "BlocksJobs.hpp"

namespace bloomrepeats {

/** Set ori so that most nucleotides have ori=1.
Alignment is inversed, if needed.
*/
class OriByMajority : public BlocksJobs {
public:
    /** Constructor */
    OriByMajority();

protected:
    void process_block_impl(Block* block, ThreadData*) const;

    const char* name_impl() const;
};

}

#endif

