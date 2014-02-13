/*
 * bloomrepeats, Find genomic repeats, using Bloom filter based prefiltration
 * Copyright (C) 2012 Boris Nagaev
 *
 * See the LICENSE file for terms of use.
 */

#ifndef BR_EXPANDER_BASE_HPP_
#define BR_EXPANDER_BASE_HPP_

#include "global.hpp"

namespace bloomrepeats {

/** Add options for fragment expansion (PairAligner) */
void add_expander_options(Processor* processor);

/** Return if two fragments can be aligned */
bool aligned(const Processor* processor,
             const Fragment& f1, const Fragment& f2);

}

#endif

