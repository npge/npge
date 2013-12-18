/*
 * bloomrepeats, Find genomic repeats, using Bloom filter based prefiltration
 * Copyright (C) 2013 Boris Nagaev
 *
 * See the LICENSE file for terms of use.
 */

#ifndef BR_BLOCK_HASH_HPP_
#define BR_BLOCK_HASH_HPP_

#include <stdint.h> // for uint32_t

#include "global.hpp"

namespace bloomrepeats {

/** Return hash of block.
Fragment::id() (i.e., sequence names, fragment positions and ori)
affects hash value.
Alignment and order of fragments does not.
*/
uint32_t block_hash(const Block* block);

}

#endif

