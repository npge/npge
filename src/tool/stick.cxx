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
#include "Filter.hpp"
#include "StickBoundaries.hpp"
#include "Output.hpp"

using namespace bloomrepeats;

class StickBoundariesPipe : public Pipe {
public:
    StickBoundariesPipe() {
        set_empty_block_set();
        add(new AddSequences);
        add(new AddBlocks);
        add(new Filter);
        add(new StickBoundaries);
        add(new Filter);
        add(new Output);
    }
};

int main(int argc, char** argv) {
    return process(argc, argv, new StickBoundariesPipe,
                   "Turn nearby fragment boundaries into one");
}

