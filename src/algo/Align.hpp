/*
 * bloomrepeats, Find genomic repeats, using Bloom filter based prefiltration
 * Copyright (C) 2012 Boris Nagaev
 *
 * See the LICENSE file for terms of use.
 */

#ifndef BR_ALIGN_HPP_
#define BR_ALIGN_HPP_

#include "Pipe.hpp"

namespace bloomrepeats {

/** Align, move and cut gaps */
class Align : public Pipe {
public:
    /** Constructor */
    Align();
};

}

#endif

