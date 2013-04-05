/*
 * bloomrepeats, Find genomic repeats, using Bloom filter based prefiltration
 * Copyright (C) 2012 Boris Nagaev
 *
 * See the LICENSE file for terms of use.
 */

#include "process.hpp"
#include "Pipe.hpp"
#include "AddBlocks.hpp"
#include "ResolveAnchors.hpp"
#include "CheckNoOverlaps.hpp"
#include "OutputPipe.hpp"

using namespace bloomrepeats;

class ResolveAnchorsPipe : public Pipe {
public:
    ResolveAnchorsPipe() {
        add(new AddBlocks, "target=other");
        add(new ResolveAnchors);
        add(new CheckNoOverlaps);
        add(new OutputPipe);
    }
};

int main(int argc, char** argv) {
    return process(argc, argv, new ResolveAnchorsPipe,
                   "Add anchors, resolve overlaps");
}

