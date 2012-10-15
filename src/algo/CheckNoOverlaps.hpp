/*
 * bloomrepeats, Find genomic repeats, using Bloom filter based prefiltration
 * Copyright (C) 2012 Boris Nagaev
 *
 * See the LICENSE file for terms of use.
 */

#ifndef BR_CHECK_NO_OVERLAPS_HPP_
#define BR_CHECK_NO_OVERLAPS_HPP_

#include "Pipe.hpp"

namespace bloomrepeats {

/** Check that there are no overlaps in block set.
BOOST_ASSERT is used.
In NDEBUG is defined, does nothing.
*/
class CheckNoOverlaps : public Pipe {
public:
    /** Constructor */
    CheckNoOverlaps();
};

}

#endif

