/*
 * NPG-explorer, Nucleotide PanGenome explorer
 * Copyright (C) 2013 Boris Nagaev
 *
 * See the LICENSE file for terms of use.
 */

#ifndef NPGE_OLD_MAKE_PRE_PANGENOME_HPP_
#define NPGE_OLD_MAKE_PRE_PANGENOME_HPP_

#include "Pipe.hpp"

namespace npge {

/** Run anchor finder, expand blocks and resolve overlaps */
class OldMakePrePangenome : public Pipe {
public:
    /** Constructor */
    OldMakePrePangenome();

protected:
    const char* name_impl() const;
};

}

#endif

