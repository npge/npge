/*
 * bloomrepeats, Find genomic repeats, using Bloom filter based prefiltration
 * Copyright (C) 2012 Boris Nagaev
 *
 * See the LICENSE file for terms of use.
 */

#ifndef BR_CONVERT_POSITION_HPP_
#define BR_CONVERT_POSITION_HPP_

#include "global.hpp"

namespace bloomrepeats {

/** Return block pos, corresponding to given fragment pos and block length */
int block_pos(const Fragment* f, int f_pos, int block_length);

}

#endif

