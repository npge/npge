/*
 * bloomrepeats, Find genomic repeats, using Bloom filter based prefiltration
 * Copyright (C) 2013 Boris Nagaev
 *
 * See the LICENSE file for terms of use.
 */

#ifndef BR_MAKE_PRE_PANGENOME_HPP_
#define BR_MAKE_PRE_PANGENOME_HPP_

#include "Pipe.hpp"

namespace bloomrepeats {

/** Run anchor finder, expand blocks and resolve overlaps */
class MakePrePangenome : public Pipe {
public:
    /** Constructor */
    MakePrePangenome();

protected:
    const char* name_impl() const;
};

}

#endif

