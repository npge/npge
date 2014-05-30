/*
 * NPG-explorer, Nucleotide PanGenome explorer
 * Copyright (C) 2013 Boris Nagaev
 *
 * See the LICENSE file for terms of use.
 */

#ifndef BR_MAKE_PANGENOME_HPP_
#define BR_MAKE_PANGENOME_HPP_

#include "Pipe.hpp"

namespace npge {

/** Run blast and Joiner until this block set becomes pangenome */
class MakePangenome : public Pipe {
public:
    /** Constructor */
    MakePangenome();

protected:
    const char* name_impl() const;
};

}

#endif

