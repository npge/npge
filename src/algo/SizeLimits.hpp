/*
 * NPG-explorer, Nucleotide PanGenome explorer
 * Copyright (C) 2012-2015 Boris Nagaev
 *
 * See the LICENSE file for terms of use.
 */

#ifndef NPGE_SIZE_LIMITS_HPP_
#define NPGE_SIZE_LIMITS_HPP_

#include "global.hpp"

namespace npge {

/** Add options for fragment and block limitations (sizes) */
void add_lite_size_limits_options(Processor* processor);

/** Add options for fragment and block limitations */
void add_size_limits_options(Processor* processor);

/** Set all parameters of fragment and block to the most permissive */
void allow_everything(Processor* processor);

}

#endif

