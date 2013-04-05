/*
 * bloomrepeats, Find genomic repeats, using Bloom filter based prefiltration
 * Copyright (C) 2012 Boris Nagaev
 *
 * See the LICENSE file for terms of use.
 */

#include "process.hpp"
#include "Pipe.hpp"
#include "AddBlocks.hpp"
#include "Connector.hpp"
#include "OverlapsResolver.hpp"
#include "CheckNoOverlaps.hpp"
#include "OutputPipe.hpp"

using namespace bloomrepeats;

class OverlapsResolverPipe : public Pipe {
public:
    OverlapsResolverPipe() {
        add(new AddBlocks);
        add(new Connector);
        add(new OverlapsResolver);
        add(new CheckNoOverlaps);
        add(new OutputPipe);
    }
};

int main(int argc, char** argv) {
    return process(argc, argv, new OverlapsResolverPipe, "Resolve overlaps");
}

