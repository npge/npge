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
#include "Connector.hpp"
#include "PrintOverlaps.hpp"

using namespace bloomrepeats;

class PrintOverlapsPipe : public Pipe {
public:
    PrintOverlapsPipe() {
        set_empty_block_set();
        add(new AddSequences);
        add(new AddBlocks);
        add(new Connector);
        add(new PrintOverlaps);
    }
};

int main(int argc, char** argv) {
    return process(argc, argv, new PrintOverlapsPipe,
                   "Print ASCII diagram with all fragments overlapping");
}

