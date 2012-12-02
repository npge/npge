/*
 * bloomrepeats, Find genomic repeats, using Bloom filter based prefiltration
 * Copyright (C) 2012 Boris Nagaev
 *
 * See the LICENSE file for terms of use.
 */

#include "process.hpp"
#include "BlockSet.hpp"
#include "Pipe.hpp"
#include "AddSequences.hpp"
#include "AddBlocks.hpp"
#include "ExternalAligner.hpp"
#include "Output.hpp"

using namespace bloomrepeats;

class ExternalAlignerPipe : public Pipe {
public:
    ExternalAlignerPipe() {
        set_block_set(boost::make_shared<BlockSet>());
        add(new AddSequences);
        add(new AddBlocks);
        add(new ExternalAligner("mafft --retree 1 --maxiterate 0 %1% > %2%"));
        add(new Output);
    }
};

int main(int argc, char** argv) {
    return process(argc, argv, new ExternalAlignerPipe,
                   "Align all blocks (by default, with fast mafft)");
}

