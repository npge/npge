/*
 * bloomrepeats, Find genomic repeats, using Bloom filter based prefiltration
 * Copyright (C) 2012 Boris Nagaev
 *
 * See the LICENSE file for terms of use.
 */

#ifndef BR_OVERLAPSLESS_UNION_HPP_
#define BR_OVERLAPSLESS_UNION_HPP_

#include "Processor.hpp"

namespace bloomrepeats {

/** Add clones of blocks from other, non-overlapping with this block set.
Added blocks are sorted
by number of fragments desc,
by length desc,
by name desc.
*/
class OverlaplessUnion : public Processor {
protected:
    bool run_impl() const;

    const char* name_impl() const;
};

}

#endif

