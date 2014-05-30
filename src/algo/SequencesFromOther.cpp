/*
 * NPG-explorer, Nucleotide PanGenome explorer
 * Copyright (C) 2014 Boris Nagaev
 *
 * See the LICENSE file for terms of use.
 */

#include "SequencesFromOther.hpp"
#include "BlockSet.hpp"

namespace npge {

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

