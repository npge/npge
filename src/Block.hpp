/*
 * bloomrepeats, Find genomic repeats, using Bloom filter based prefiltration
 * Copyright (C) 2012 Boris Nagaev
 *
 * See the LICENSE file for terms of use.
 */

#ifndef BR_BLOCK_HPP_
#define BR_BLOCK_HPP_

#include <set>

#include "global.hpp"

namespace bloomrepeats {

class Block : public std::set<FragmentPtr> {
};

}

#endif

