/*
 * NPG-explorer, Nucleotide PanGenome explorer
 * Copyright (C) 2012 Boris Nagaev
 *
 * See the LICENSE file for terms of use.
 */

#ifndef BR_SIZE_LIMITS_HPP_
#define BR_SIZE_LIMITS_HPP_

#include "global.hpp"

namespace npge {

/** Add options for fragment and block limitations */
void add_size_limits_options(Processor* processor);

/** Set all parameters of fragment and block to the most permissive */
void allow_everything(Processor* processor);

}

#endif

