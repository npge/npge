/*
 * NPG-explorer, Nucleotide PanGenome explorer
 * Copyright (C) 2013 Boris Nagaev
 *
 * See the LICENSE file for terms of use.
 */

#ifndef NPGE_MKDIR_HPP_
#define NPGE_MKDIR_HPP_

#include "Processor.hpp"

namespace npge {

/** Create directory */
class MkDir : public Processor {
public:
    /** Constructor */
    MkDir();

protected:
    void run_impl() const;
    const char* name_impl() const;
};

}

#endif

