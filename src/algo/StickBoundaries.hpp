/*
 * bloomrepeats, Find genomic repeats, using Bloom filter based prefiltration
 * Copyright (C) 2012 Boris Nagaev
 *
 * See the LICENSE file for terms of use.
 */

#ifndef BR_STICK_BOUNDARIES_HPP_
#define BR_STICK_BOUNDARIES_HPP_

#include "Processor.hpp"
#include "config.hpp"

namespace bloomrepeats {

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
    StickBoundaries(int min_distance = BOUNDARIES_MIN_DISTANCE);

protected:
    void run_impl() const;

    const char* name_impl() const;
};

}

#endif

