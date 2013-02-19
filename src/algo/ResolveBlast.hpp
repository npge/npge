/*
 * bloomrepeats, Find genomic repeats, using Bloom filter based prefiltration
 * Copyright (C) 2012 Boris Nagaev
 *
 * See the LICENSE file for terms of use.
 */

#ifndef BR_RESOLVE_BLAST_HPP_
#define BR_RESOLVE_BLAST_HPP_

#include "Pipe.hpp"

namespace bloomrepeats {

/** AddRest, ConSeq, AddBlastBlocks, OverlapsResolver2, DeConSeq.
Source block set is used as input.
*/
class ResolveBlast : public Pipe {
public:
    /** Constructor */
    ResolveBlast(BlockSetPtr source = BlockSetPtr());
};

}

#endif

