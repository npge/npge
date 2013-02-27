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
#include "Connector.hpp"
#include "Rest.hpp"
#include "CheckNoOverlaps.hpp"
#include "OutputPipe.hpp"

using namespace bloomrepeats;

class RestPipe : public Pipe {
public:
    RestPipe() {
        set_empty_block_set();
        set_empty_other();
        add(new AddSequences, "target=other");
        add(new AddBlocks, "target=other");
        add(new Connector, "target=other");
        add(new Rest);
        add(new CheckNoOverlaps);
        add(new OutputPipe);
    }
};

int main(int argc, char** argv) {
    return process(argc, argv, new RestPipe, "Select unique fragments");
}

