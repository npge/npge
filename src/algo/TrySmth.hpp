/*
 * bloomrepeats, Find genomic repeats, using Bloom filter based prefiltration
 * Copyright (C) 2014 Boris Nagaev
 *
 * See the LICENSE file for terms of use.
 */

#ifndef BR_TRY_SMTH_HPP
#define BR_TRY_SMTH_HPP

#include "Pipe.hpp"

namespace bloomrepeats {

/** Try to do something, align (and filter), restore original if worse */
class TrySmth : public Pipe {
public:
    /** Constructor */
    TrySmth();

protected:
    const char* name_impl() const;
};

}

#endif

