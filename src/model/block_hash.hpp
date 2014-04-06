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

/** Return hash of block.
Hashes of blocks of blockset are XOR'ed.
Blocks of <=1 fragment are skipped.
*/
uint32_t blockset_hash(const BlockSet& block_set, int workers = 1);

/** Return block id (<size>x<length>) */
std::string block_id(const Block* block);

/** Return if block has fragments from same genome */
bool has_repeats(const Block* block);

/** Return number of genomes occupied by block set */
int genomes_number(const BlockSet& block_set);

/** Return block name.
Block name format:
 - first letter
   - r for repeat (>= 2 fragments on same genome)
   - s for exact steam (= 1 fragment on each genome)
   - h for other blocks of >=2 fragments
   - u for blocks of 1 fragment.
 - block_id
*/
std::string block_name(const Block* block, int genomes);

}

#endif

