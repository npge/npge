/*
 * bloomrepeats, Find genomic repeats, using Bloom filter based prefiltration
 * Copyright (C) 2012 Boris Nagaev
 *
 * See the LICENSE file for terms of use.
 */

#include "SequencesFromOther.hpp"
#include "BlockSet.hpp"

namespace bloomrepeats {

SequencesFromOther::SequencesFromOther(const BlockSetPtr& source) {
    set_other(source);
    declare_bs("other", "Source");
    declare_bs("target", "Destination");
}

void SequencesFromOther::run_impl() const {
    block_set()->add_sequences(other()->seqs());
}

const char* SequencesFromOther::name_impl() const {
    return "Copy sequences from other block set";
}

}

