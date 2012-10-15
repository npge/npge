/*
 * bloomrepeats, Find genomic repeats, using Bloom filter based prefiltration
 * Copyright (C) 2012 Boris Nagaev
 *
 * See the LICENSE file for terms of use.
 */

#include "Swap.hpp"
#include "BlockSet.hpp"

namespace bloomrepeats {

Swap::Swap(const BlockSetPtr& other):
    OtherBlockSet(other)
{ }

bool Swap::run_impl() const {
    block_set()->swap(*other());
    return true;
}

}

