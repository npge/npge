/*
 * bloomrepeats, Find genomic repeats, using Bloom filter based prefiltration
 * Copyright (C) 2012 Boris Nagaev
 *
 * See the LICENSE file for terms of use.
 */

#ifndef BR_UNIQUE_NAMES_HPP_
#define BR_UNIQUE_NAMES_HPP_

#include "Processor.hpp"

namespace bloomrepeats {

/** Set unique names to all blocks of this block set.
Firstly, Block::set_name_from_fragments() is used, if name is null.
Then Block::set_random_name() is called untill the name is unique.
*/
class UniqueNames : public Processor {
protected:
    /** Apply the action */
    bool run_impl() const;
};

}

#endif

