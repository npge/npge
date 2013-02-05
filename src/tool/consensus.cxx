/*
 * bloomrepeats, Find genomic repeats, using Bloom filter based prefiltration
 * Copyright (C) 2012 Boris Nagaev
 *
 * See the LICENSE file for terms of use.
 */

#include "process.hpp"
#include "Pipe.hpp"
#include "AddSequences.hpp"
#include "AddBlocks.hpp"
#include "ExternalAligner.hpp"
#include "Consensus.hpp"

using namespace bloomrepeats;

class ExternalAlignerPipe : public Pipe {
public:
    ExternalAlignerPipe() {
        set_empty_block_set();
        add(new AddSequences);
        add(new AddBlocks);
        add(new ExternalAligner);
        Consensus* consensus = new Consensus;
        consensus->set_remove_after(false);
        add(consensus);
    }
};

int main(int argc, char** argv) {
    return process(argc, argv, new ExternalAlignerPipe,
                   "Align all blocks and write consensuses to file");
}

