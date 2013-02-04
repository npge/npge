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
#include "UniqueNames.hpp"
#include "CheckNoOverlaps.hpp"
#include "Output.hpp"

using namespace bloomrepeats;

class RestPipe : public Pipe {
public:
    RestPipe() {
        set_empty_block_set();
        add(new AddSequences);
        add(new AddBlocks);
        add(new Connector);
        add(new Rest, THIS_TO_OTHER | THIS_TO_THIS);
        add(new CheckNoOverlaps);
        add(new UniqueNames);
        add(new Output);
    }
};

int main(int argc, char** argv) {
    return process(argc, argv, new RestPipe,
                   "Add unique fragments to current block set");
}

