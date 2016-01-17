/*
 * NPG-explorer, Nucleotide PanGenome explorer
 * Copyright (C) 2012-2016 Boris Nagaev
 *
 * See the LICENSE file for terms of use.
 */

#ifndef NPGE_REFINE_ALIGNMENT_HPP_
#define NPGE_REFINE_ALIGNMENT_HPP_

#include "global.hpp"

namespace npge {

/** Move mismatches and remove pure gap columns */
void refine_alignment(Strings& aligned);

}

#endif

