/*
 * NPG-explorer, Nucleotide PanGenome explorer
 * Copyright (C) 2012 Boris Nagaev
 *
 * See the LICENSE file for terms of use.
 */

#ifndef NPGE_EXPANDER_BASE_HPP_
#define NPGE_EXPANDER_BASE_HPP_

#include "global.hpp"

namespace npge {

/** Add options for fragment expansion (PairAligner) */
void add_expander_options(Processor* processor);

/** Apply options to pair aligner */
void apply_pair_aligner_options(PairAligner* pa, const Processor* p);

/** Return if two fragments can be aligned */
bool aligned(const Processor* processor,
             const Fragment& f1, const Fragment& f2);

}

#endif

