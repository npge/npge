/*
 * NPG-explorer, Nucleotide PanGenome explorer
 * Copyright (C) 2014 Boris Nagaev
 *
 * See the LICENSE file for terms of use.
 */

#ifndef NPGE_CLEAN_UP_HPP_
#define NPGE_CLEAN_UP_HPP_

#include "Pipe.hpp"

namespace npge {

/** Connect, resolve overlaps, expand, filter (FIXME) */
class CleanUp : public Pipe {
public:
    /** Constructor */
    CleanUp();

protected:
    const char* name_impl() const;
};

}

#endif

