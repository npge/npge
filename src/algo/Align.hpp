/*
 * NPG-explorer, Nucleotide PanGenome explorer
 * Copyright (C) 2014 Boris Nagaev
 *
 * See the LICENSE file for terms of use.
 */

#ifndef NPGE_ALIGN_HPP_
#define NPGE_ALIGN_HPP_

#include "Pipe.hpp"

namespace npge {

/** Align, move and cut gaps */
class LiteAlign : public Pipe {
public:
    /** Constructor */
    LiteAlign();

protected:
    const char* name_impl() const;
};

/** LiteAlign + Filter + SelfOverlapsResolver */
class Align : public Pipe {
public:
    /** Constructor */
    Align();

protected:
    const char* name_impl() const;
};

}

#endif

