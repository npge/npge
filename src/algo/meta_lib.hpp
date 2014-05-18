/*
 * bloomrepeats, Find genomic repeats, using Bloom filter based prefiltration
 * Copyright (C) 2014 Boris Nagaev
 *
 * See the LICENSE file for terms of use.
 */

#ifndef BR_META_LIB_HPP_
#define BR_META_LIB_HPP_

#include "global.hpp"

namespace bloomrepeats {

/** Add stabdard processors to Meta.
This function is called from constructor of Meta.
*/
void add_meta_lib(Meta* meta);

}

#endif

