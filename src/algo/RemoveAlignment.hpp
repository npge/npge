/*
 * NPG-explorer, Nucleotide PanGenome explorer
 * Copyright (C) 2014 Boris Nagaev
 *
 * See the LICENSE file for terms of use.
 */

#ifndef NPGE_REMOVE_ALIGNMENT_HPP_
#define NPGE_REMOVE_ALIGNMENT_HPP_

#include "Processor.hpp"

namespace npge {

/** Remove alignment rows of all fragments */
class RemoveAlignment : public Processor {
public:
    /** Constructor */
    RemoveAlignment();

protected:
    void run_impl() const;
    const char* name_impl() const;
};

}

#endif

