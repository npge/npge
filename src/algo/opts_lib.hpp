/*
 * bloomrepeats, Find genomic repeats, using Bloom filter based prefiltration
 * Copyright (C) 2014 Boris Nagaev
 *
 * See the LICENSE file for terms of use.
 */

#ifndef BR_OPTS_LIB_HPP_
#define BR_OPTS_LIB_HPP_

#include "global.hpp"

namespace bloomrepeats {

/** Add global options to Meta */
void add_opts(Meta* meta);

}

#endif

