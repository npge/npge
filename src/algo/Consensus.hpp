/*
 * NPG-explorer, Nucleotide PanGenome explorer
 * Copyright (C) 2014 Boris Nagaev
 *
 * See the LICENSE file for terms of use.
 */

#ifndef NPGE_CONSENSUS_HPP_
#define NPGE_CONSENSUS_HPP_

#include "Pipe.hpp"

namespace npge {

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

