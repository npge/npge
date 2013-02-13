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
}

bool SequencesFromOther::run_impl() const {
    block_set()->add_sequences(other()->seqs());
}

}

