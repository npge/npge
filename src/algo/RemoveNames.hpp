/*
 * bloomrepeats, Find genomic repeats, using Bloom filter based prefiltration
 * Copyright (C) 2013 Boris Nagaev
 *
 * See the LICENSE file for terms of use.
 */

#ifndef BR_REMOVE_NAMES_HPP_
#define BR_REMOVE_NAMES_HPP_

#include "Processor.hpp"

namespace bloomrepeats {

/** Set names of all sequences/blocks to empty */
class RemoveNames : public Processor {
public:
    /** Constructor */
    RemoveNames();

protected:
    void run_impl() const;

    const char* name_impl() const;
};

}

#endif

