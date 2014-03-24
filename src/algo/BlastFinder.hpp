/*
 * bloomrepeats, Find genomic repeats, using Bloom filter based prefiltration
 * Copyright (C) 2012 Boris Nagaev
 *
 * See the LICENSE file for terms of use.
 */

#ifndef BR_BLAST_FINDER_HPP_
#define BR_BLAST_FINDER_HPP_

#include "Pipe.hpp"

namespace bloomrepeats {

/** Takes sequences, add blast hits as blocks of these sequences.
Run BlastFinder on sequence-only blockset (without blocks).
It adds blast hits as blocks to this blockset.
*/
class BlastFinder : public Pipe {
public:
    /** Constructor */
    BlastFinder();

private:
    mutable std::string consensus_;
    mutable std::string hits_;

    std::string consensus() const;
    std::string hits() const;
};

}

#endif

