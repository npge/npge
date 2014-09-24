/*
 * NPG-explorer, Nucleotide PanGenome explorer
 * Copyright (C) 2014 Boris Nagaev
 *
 * See the LICENSE file for terms of use.
 */

#ifndef NPGE_OVERLAPS_RESOLVER2_HPP_
#define NPGE_OVERLAPS_RESOLVER2_HPP_

#include "OverlapsResolver.hpp"

namespace npge {

/** Resolve overlaping fragments (version 2).
\warning Source blocks are taken from "other" blockset.
    Other blockset can be changed (similar to StickBoundaries).
*/
class OverlapsResolver2 : public OverlapsResolver {
public:
    /** Constructor */
    OverlapsResolver2();

protected:
    void run_impl() const;

    const char* name_impl() const;
};

}

#endif

