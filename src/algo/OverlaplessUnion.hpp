/*
 * NPG-explorer, Nucleotide PanGenome explorer
 * Copyright (C) 2012-2015 Boris Nagaev
 *
 * See the LICENSE file for terms of use.
 */

#ifndef NPGE_OVERLAPSLESS_UNION_HPP_
#define NPGE_OVERLAPSLESS_UNION_HPP_

#include "Processor.hpp"

namespace npge {

/** Add clones of blocks from other, non-overlapping with this blockset.
Added blocks are sorted
by number of fragments desc,
by length desc,
by name desc.
*/
class OverlaplessUnion : public Processor {
public:
    /** Constructor */
    OverlaplessUnion();

protected:
    void run_impl() const;

    const char* name_impl() const;
};

}

#endif

