/*
 * NPG-explorer, Nucleotide PanGenome explorer
 * Copyright (C) 2012 Boris Nagaev
 *
 * See the LICENSE file for terms of use.
 */

#ifndef NPGE_CLEAR_HPP_
#define NPGE_CLEAR_HPP_

#include "Processor.hpp"

namespace npge {

/** Remove all blocks and/or sequences from blockset */
class Clear : public Processor {
public:
    /** Constructor */
    Clear();

protected:
    void run_impl() const;

    const char* name_impl() const;
};

}

#endif

