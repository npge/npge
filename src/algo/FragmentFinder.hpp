/*
 * NPG-explorer, Nucleotide PanGenome explorer
 * Copyright (C) 2012-2016 Boris Nagaev
 *
 * See the LICENSE file for terms of use.
 */

#ifndef NPGE_FRAGMENT_FINDER_HPP_
#define NPGE_FRAGMENT_FINDER_HPP_

#include "global.hpp"
#include "Processor.hpp"

namespace npge {

/** Locate fragment by its sequence */
class FragmentFinder : public Processor {
public:
    /** Default constructor */
    FragmentFinder();

protected:
    void run_impl() const;

    const char* name_impl() const;
};

}

#endif

