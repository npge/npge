/*
 * bloomrepeats, Find genomic repeats, using Bloom filter based prefiltration
 * Copyright (C) 2012 Boris Nagaev
 *
 * See the LICENSE file for terms of use.
 */

#include "process.hpp"
#include "BlockSet.hpp"
#include "Pipe.hpp"
#include "AddBlocks.hpp"
#include "Filter.hpp"
#include "StickBoundaries.hpp"
#include "OutputPipe.hpp"

using namespace bloomrepeats;

class StickBoundariesPipe : public Pipe {
public:
    StickBoundariesPipe() {
        add(new AddBlocks);
        add(new Filter);
        add(new StickBoundaries);
        add(new Filter);
        add(new OutputPipe);
    }
};

int main(int argc, char** argv) {
    return process(argc, argv, new StickBoundariesPipe,
                   "Turn nearby fragment boundaries into one");
}

