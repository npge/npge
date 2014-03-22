/*
 * bloomrepeats, Find genomic repeats, using Bloom filter based prefiltration
 * Copyright (C) 2012 Boris Nagaev
 *
 * See the LICENSE file for terms of use.
 */

#ifndef BR_CLEAR_HPP_
#define BR_CLEAR_HPP_

#include "Processor.hpp"

namespace bloomrepeats {

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

