/*
 * bloomrepeats, Find genomic repeats, using Bloom filter based prefiltration
 * Copyright (C) 2014 Boris Nagaev
 *
 * See the LICENSE file for terms of use.
 */

#ifndef BR_REFINE_ALIGNMENT_HPP_
#define BR_REFINE_ALIGNMENT_HPP_

#include "global.hpp"

namespace bloomrepeats {

/** Move mismatches and remove pure gap columns */
void refine_alignment(Strings& aligned);

}

#endif

