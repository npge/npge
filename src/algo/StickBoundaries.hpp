/*
 * NPG-explorer, Nucleotide PanGenome explorer
 * Copyright (C) 2012 Boris Nagaev
 *
 * See the LICENSE file for terms of use.
 */

#ifndef NPGE_STICK_BOUNDARIES_HPP_
#define NPGE_STICK_BOUNDARIES_HPP_

#include "Processor.hpp"

namespace npge {

/** Move boundaries of fragments to turn nearby boundaries into one.
It is recommended to \ref Filter "filter" fragments before
applying this processor.

\warning Connections between fragments may be broken.
    They must be \ref Connector "re-connected" after applying this processor.

\warning Alignment of blocks may be broken.
*/
class StickBoundaries : public Processor {
public:
    /** Constructor */
    StickBoundaries();

protected:
    void run_impl() const;

    const char* name_impl() const;
};

}

#endif

