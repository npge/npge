/*
 * NPG-explorer, Nucleotide PanGenome explorer
 * Copyright (C) 2014 Boris Nagaev
 *
 * See the LICENSE file for terms of use.
 */

#ifndef NPGE_BLOCK_FINDER_HPP_
#define NPGE_BLOCK_FINDER_HPP_

#include "global.hpp"
#include "Processor.hpp"

namespace npge {

class FragmentFinder;
class OverlapFinder;

/** Locate block by sequence of its fragments */
class BlockFinder : public Processor {
public:
    /** Default constructor */
    BlockFinder();

protected:
    void run_impl() const;

    const char* name_impl() const;

private:
    FragmentFinder* ff_;
    OverlapFinder* of_;
};

}

#endif

