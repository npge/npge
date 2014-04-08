/*
 * bloomrepeats, Find genomic repeats, using Bloom filter based prefiltration
 * Copyright (C) 2012 Boris Nagaev
 *
 * See the LICENSE file for terms of use.
 */

#ifndef BR_CLEAN_UP_HPP_
#define BR_CLEAN_UP_HPP_

#include "Pipe.hpp"

namespace bloomrepeats {

/** Connect, resolve overlaps, expand, filter (FIXME) */
class CleanUp : public Pipe {
public:
    /** Constructor */
    CleanUp();

protected:
    const char* name_impl() const;
};

}

#endif

