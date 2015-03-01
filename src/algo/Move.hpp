/*
 * NPG-explorer, Nucleotide PanGenome explorer
 * Copyright (C) 2012-2015 Boris Nagaev
 *
 * See the LICENSE file for terms of use.
 */

#ifndef NPGE_MOVE_HPP_
#define NPGE_MOVE_HPP_

#include "Processor.hpp"

namespace npge {

/** Move all blocks from other blockset to target blockset */
class Move : public Processor {
public:
    /** Constructor */
    Move();

protected:
    void run_impl() const;
    const char* name_impl() const;
};

}

#endif

