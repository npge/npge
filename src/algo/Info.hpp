/*
 * NPG-explorer, Nucleotide PanGenome explorer
 * Copyright (C) 2012 Boris Nagaev
 *
 * See the LICENSE file for terms of use.
 */

#ifndef NPGE_INFO_HPP_
#define NPGE_INFO_HPP_

#include "Processor.hpp"

namespace npge {

class Stats;

/** Print human readable summary and statistics about block set (+ subsets) */
class Info : public Processor {
public:
    /** Constructor */
    Info();

protected:
    void run_impl() const;

    const char* name_impl() const;

private:
    Stats* stats_;
};

}

#endif

