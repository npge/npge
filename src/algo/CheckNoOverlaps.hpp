/*
 * NPG-explorer, Nucleotide PanGenome explorer
 * Copyright (C) 2012 Boris Nagaev
 *
 * See the LICENSE file for terms of use.
 */

#ifndef NPGE_CHECK_NO_OVERLAPS_HPP_
#define NPGE_CHECK_NO_OVERLAPS_HPP_

#include "Pipe.hpp"

namespace npge {

/** Check that there are no overlaps in block set.
Throws an exception if overlaps are found.
*/
class CheckNoOverlaps : public Pipe {
public:
    /** Constructor */
    CheckNoOverlaps();

protected:
    const char* name_impl() const;
};

}

#endif

