/*
 * NPG-explorer, Nucleotide PanGenome explorer
 * Copyright (C) 2014 Boris Nagaev
 *
 * See the LICENSE file for terms of use.
 */

#ifndef BR_FIX_ENDS_HPP_
#define BR_FIX_ENDS_HPP_

#include "BlocksJobs.hpp"

namespace npge {

/** Cut bad aligned ends.

Alignment is preserved.

Slice subblock such that it's both ends of length min_fragment
are good. Good column is an ident column without gaps.
Percentage of good columns is greater or equal to
min_identity
Start colum is good.
It tries to find exact boundary of good subblock by
selecting better boundary in min_fragment * (1 - min_identity)
columns after first good boundary.

Remove blocks failed to fix.
*/
class FixEnds : public BlocksJobs {
public:
    /** Constructor */
    FixEnds();

protected:
    ThreadData* before_thread_impl() const;
    void process_block_impl(Block* block, ThreadData*) const;
    void after_thread_impl(ThreadData* data) const;
    const char* name_impl() const;
};

}

#endif

