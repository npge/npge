/*
 * NPG-explorer, Nucleotide PanGenome explorer
 * Copyright (C) 2014 Boris Nagaev
 *
 * See the LICENSE file for terms of use.
 */

#ifndef BR_EXACT_STEM_BSA_HPP_
#define BR_EXACT_STEM_BSA_HPP_

#include "Processor.hpp"

namespace npge {

/** Replace all non-stem blocks with gaps in block set alignment */
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

