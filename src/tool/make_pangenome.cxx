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
#include "CleanUp.hpp"
#include "CheckNoOverlaps.hpp"
#include "OutputPipe.hpp"

using namespace bloomrepeats;

class CleanUpPipe : public Pipe {
public:
    CleanUpPipe() {
        add(new AddSequences);
        add(new AddBlocks);
        add(new Connector);
        add(new CleanUp);
        add(new CheckNoOverlaps);
        add(new OutputPipe);
    }
};

int main(int argc, char** argv) {
    return process(argc, argv, new CleanUpPipe,
                   "Resolve overlaps and expand blocks");
}

