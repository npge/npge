/*
 * bloomrepeats, Find genomic repeats, using Bloom filter based prefiltration
 * Copyright (C) 2012 Boris Nagaev
 *
 * See the LICENSE file for terms of use.
 */

#include "AddBlastBlocks.hpp"
#include "SequencesFromOther.hpp"
#include "Clear.hpp"
#include "ConSeq.hpp"
#include "BlastFinder.hpp"
#include "DeConSeq.hpp"

namespace bloomrepeats {

AddBlastBlocks::AddBlastBlocks(BlockSetPtr source):
    Pipe(source) {
    add(new SequencesFromOther);
    add(new Clear, "target=consensus --clear-seqs:=1");
    add(new ConSeq, "target=consensus other=other");
    add(new BlastFinder, "target=consensus");
    add(new DeConSeq, "target=target other=consensus");
}

}

