/*
 * NPG-explorer, Nucleotide PanGenome explorer
 * Copyright (C) 2012-2015 Boris Nagaev
 *
 * See the LICENSE file for terms of use.
 */

#ifndef NPGE_REST_HPP_
#define NPGE_REST_HPP_

#include "Processor.hpp"

namespace npge {

/** Return new blockset of blocks of nucleotides, not included in this set.
From sequences involved in this blockset, nucleotides are selected,
not included in this blockset. They are grouped into fragments.
Each fragment is inserted into one block.
These blocks are inserted into resulting blockset.
*/
class Rest : public Processor {
public:
    /** Constructor */
    Rest(const BlockSetPtr& source = BlockSetPtr());

protected:
    void run_impl() const;
    const char* name_impl() const;
};

}

#endif

