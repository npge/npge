/*
 * NPG-explorer, Nucleotide PanGenome explorer
 * Copyright (C) 2013 Boris Nagaev
 *
 * See the LICENSE file for terms of use.
 */

#ifndef NPGE_OLD_MAKE_PANGENOME_HPP_
#define NPGE_OLD_MAKE_PANGENOME_HPP_

#include "Pipe.hpp"

namespace npge {

/** Run blast and Joiner until this block set becomes pangenome */
class OldMakePangenome : public Pipe {
public:
    /** Constructor */
    OldMakePangenome();

protected:
    const char* name_impl() const;
};

}

#endif
