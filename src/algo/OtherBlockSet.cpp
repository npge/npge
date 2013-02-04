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
    other_(other), processor_(0), other_block_set_(0)
{ }

BlockSetPtr OtherBlockSet::other() const {
    if (other_) {
        return other_;
    }
    if (processor_) {
        return processor_->block_set();
    }
    if (other_block_set_) {
        return other_block_set_->other();
    }
    return BlockSetPtr();
}

void OtherBlockSet::set_other(const BlockSetPtr& other) {
    other_ = other;
    processor_ = 0;
    other_block_set_ = 0;
}

void OtherBlockSet::set_processor(const Processor* processor) {
    other_.reset();
    processor_ = processor;
    other_block_set_ = 0;
}

void OtherBlockSet::set_other_block_set(OtherBlockSet* other_block_set) {
    other_.reset();
    processor_ = 0;
    other_block_set_ = other_block_set;
}

}

