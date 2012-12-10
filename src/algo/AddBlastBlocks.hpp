/*
 * bloomrepeats, Find genomic repeats, using Bloom filter based prefiltration
 * Copyright (C) 2012 Boris Nagaev
 *
 * See the LICENSE file for terms of use.
 */

#ifndef BR_ADD_BLAST_BLOCKS_HPP_
#define BR_ADD_BLAST_BLOCKS_HPP_

#include "Pipe.hpp"

namespace bloomrepeats {

/** Make consensuses, run blast and read blast hits.
It is recommended to run CleanUp after this processor.
*/
class AddBlastBlocks : public Pipe {
public:
    /** Constructor */
    AddBlastBlocks();
};

}

#endif

