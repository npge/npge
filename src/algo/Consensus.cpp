/*
 * NPG-explorer, Nucleotide PanGenome explorer
 * Copyright (C) 2014 Boris Nagaev
 *
 * See the LICENSE file for terms of use.
 */

#include "Consensus.hpp"
#include "ConSeq.hpp"
#include "Output.hpp"
#include "Clear.hpp"

namespace npge {

Consensus::Consensus() {
    add(new ConSeq, "target=cons other=target");
    add(new Output, "target=cons prefix|cons- "
        "--cons-dump-seq:=1 --cons-dump-block:=0");
    add(new Clear, "target=cons --clear-blocks:=1 --clear-seqs:=1");
    declare_bs("target", "Target blockset");
}

const char* Consensus::name_impl() const {
    return "Consensus writer";
}

}

