/*
 * NPG-explorer, Nucleotide PanGenome explorer
 * Copyright (C) 2014 Boris Nagaev
 *
 * See the LICENSE file for terms of use.
 */

#ifndef NPGE_BLOCK_SET_ALIGNMENT_PROCESSOR_HPP_
#define NPGE_BLOCK_SET_ALIGNMENT_PROCESSOR_HPP_

#include "global.hpp"
#include "Processor.hpp"

namespace npge {

/** Build blockset alignment */
class FindBSA : public Processor {
public:
    /** Constructor */
    FindBSA();

protected:
    void run_impl() const;

    const char* name_impl() const;
};

}

#endif

