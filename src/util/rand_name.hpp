/*
 * NPG-explorer, Nucleotide PanGenome explorer
 * Copyright (C) 2013 Boris Nagaev
 *
 * See the LICENSE file for terms of use.
 */

#ifndef NPGE_RAND_NAME_HPP_
#define NPGE_RAND_NAME_HPP_

#include <string>

namespace npge {

/** Return random seed */
int make_seed();

/** Return random name */
std::string rand_name(int size);

}

#endif

