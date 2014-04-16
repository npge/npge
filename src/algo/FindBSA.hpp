/*
 * bloomrepeats, Find genomic repeats, using Bloom filter based prefiltration
 * Copyright (C) 2014 Boris Nagaev
 *
 * See the LICENSE file for terms of use.
 */

#ifndef BR_BLOCK_SET_ALIGNMENT_PROCESSOR_HPP_
#define BR_BLOCK_SET_ALIGNMENT_PROCESSOR_HPP_

#include "global.hpp"
#include "Processor.hpp"

namespace bloomrepeats {

/** Build block set alignment */
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

