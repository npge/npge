/*
 * NPG-explorer, Nucleotide PanGenome explorer
 * Copyright (C) 2014 Boris Nagaev
 *
 * See the LICENSE file for terms of use.
 */

#ifndef NPGE_MERGE_UNIQUE_HPP_
#define NPGE_MERGE_UNIQUE_HPP_

#include "Processor.hpp"

namespace npge {

/** Merge unique fragments with common neighbours into blocks */
class MergeUnique : public Processor {
public:
    /** Constructor */
    MergeUnique();

protected:
    void run_impl() const;

    const char* name_impl() const;
};

}

#endif

