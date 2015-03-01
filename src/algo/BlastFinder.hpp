/*
 * NPG-explorer, Nucleotide PanGenome explorer
 * Copyright (C) 2012-2015 Boris Nagaev
 *
 * See the LICENSE file for terms of use.
 */

#ifndef NPGE_BLAST_FINDER_HPP_
#define NPGE_BLAST_FINDER_HPP_

#include "Pipe.hpp"

namespace npge {

/** Takes sequences, add blast hits as blocks of these sequences.
Run BlastFinder on sequence-only blockset (without blocks).
It adds blast hits as blocks to this blockset.
*/
class BlastFinder : public Pipe {
public:
    /** Constructor */
    BlastFinder();

protected:
    void run_impl() const;

private:
    mutable std::string consensus_;
    mutable std::string hits_;

    std::string consensus() const;
    std::string hits() const;
};

}

#endif

