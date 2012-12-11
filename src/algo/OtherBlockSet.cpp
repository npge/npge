/*
 * bloomrepeats, Find genomic repeats, using Bloom filter based prefiltration
 * Copyright (C) 2012 Boris Nagaev
 *
 * See the LICENSE file for terms of use.
 */

#include "OtherBlockSet.hpp"
#include "Processor.hpp"

namespace bloomrepeats {

OtherBlockSet::OtherBlockSet(const BlockSetPtr& other):
    other_(other)
{ }

OtherBlockSet::OtherBlockSet(const Processor* processor):
    processor_(processor)
{ }

BlockSetPtr OtherBlockSet::other() const {
    return other_ ? : processor_->block_set() ? : BlockSetPtr();
}

}

