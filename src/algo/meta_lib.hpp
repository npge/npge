/*
 * NPG-explorer, Nucleotide PanGenome explorer
 * Copyright (C) 2012-2015 Boris Nagaev
 *
 * See the LICENSE file for terms of use.
 */

#ifndef NPGE_META_LIB_HPP_
#define NPGE_META_LIB_HPP_

#include "global.hpp"

namespace npge {

/** Add stabdard processors to Meta.
This function is called from constructor of Meta.
*/
void add_meta_lib(Meta* meta);

}

#endif

