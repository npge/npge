/*
 * bloomrepeats, Find genomic repeats, using Bloom filter based prefiltration
 * Copyright (C) 2012 Boris Nagaev
 *
 * See the LICENSE file for terms of use.
 */

#ifndef BR_INFO_HPP_
#define BR_INFO_HPP_

#include "Processor.hpp"

namespace bloomrepeats {

/** Print human readable summary and statistics about block set (+ subsets) */
class Info : public Processor {
public:
    /** Constructor */
    Info();

protected:
    bool run_impl() const;

    const char* name_impl() const;

private:
    Stats* stats_;
};

}

#endif

