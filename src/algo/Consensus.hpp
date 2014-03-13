/*
 * bloomrepeats, Find genomic repeats, using Bloom filter based prefiltration
 * Copyright (C) 2012 Boris Nagaev
 *
 * See the LICENSE file for terms of use.
 */

#ifndef BR_CONSENSUS_HPP_
#define BR_CONSENSUS_HPP_

#include "Pipe.hpp"

namespace bloomrepeats {

/** Write consensuses of all blocks to a file */
class Consensus : public Pipe {
public:
    /** Constructor */
    Consensus();

protected:
    const char* name_impl() const;
};

}

#endif

