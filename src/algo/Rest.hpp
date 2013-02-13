/*
 * bloomrepeats, Find genomic repeats, using Bloom filter based prefiltration
 * Copyright (C) 2012 Boris Nagaev
 *
 * See the LICENSE file for terms of use.
 */

#ifndef BR_REST_HPP_
#define BR_REST_HPP_

#include "Processor.hpp"

namespace bloomrepeats {

/** Return new block set of blocks of nucleotides, not included in this set.
From sequences involved in this block set, nucleotides are selected,
not included in this block set. They are grouped into fragments.
Each fragment is inserted into one block.
These blocks are inserted into resulting block set.
\warning Fragments must be \ref Connector "connected"
   for this to work correctly.
*/
class Rest : public Processor {
public:
    /** Constructor */
    Rest(const BlockSetPtr& source = BlockSetPtr());

protected:
    /** Apply the action */
    bool run_impl() const;
};

}

#endif

