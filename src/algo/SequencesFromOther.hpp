/*
 * bloomrepeats, Find genomic repeats, using Bloom filter based prefiltration
 * Copyright (C) 2012 Boris Nagaev
 *
 * See the LICENSE file for terms of use.
 */

#ifndef BR_SEQUENCES_FROM_OTHER_HPP_
#define BR_SEQUENCES_FROM_OTHER_HPP_

#include "Processor.hpp"

namespace bloomrepeats {

/** Copy sequences from other block set.
Wrapper of BlockSet::add_sequences().
*/
class SequencesFromOther : public Processor {
public:
    /** Constructor
    \param source BlockSet, from which sequences will be added to block_set().
    */
    SequencesFromOther(const BlockSetPtr& source = BlockSetPtr());

protected:
    /** Apply the action */
    void run_impl() const;
};

}

#endif

