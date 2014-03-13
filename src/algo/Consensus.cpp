/*
 * bloomrepeats, Find genomic repeats, using Bloom filter based prefiltration
 * Copyright (C) 2012 Boris Nagaev
 *
 * See the LICENSE file for terms of use.
 */

#include "Consensus.hpp"
#include "ConSeq.hpp"
#include "Output.hpp"
#include "Clear.hpp"

namespace bloomrepeats {

Consensus::Consensus() {
    add(new ConSeq, "target=cons other=target");
    add(new Output, "target=cons prefix|cons- "
        "--cons-dump-seq:=1 --cons-dump-block:=0");
    add(new Clear, "target=cons --clear-blocks:=1 --clear-seqs:=1");
}

const char* Consensus::name_impl() const {
    return "Consensus writer";
}

}

