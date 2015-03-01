/*
 * NPG-explorer, Nucleotide PanGenome explorer
 * Copyright (C) 2012-2015 Boris Nagaev
 *
 * See the LICENSE file for terms of use.
 */

#ifndef NPGE_EXACT_STEM_BSA_HPP_
#define NPGE_EXACT_STEM_BSA_HPP_

#include "Processor.hpp"

namespace npge {

/** Replace all non-stem blocks with gaps in blockset alignment */
class ExactStemBSA : public Processor {
public:
    /** Constructor */
    ExactStemBSA();

protected:
    void run_impl() const;

    const char* name_impl() const;
};

}

#endif

