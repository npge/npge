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
#include "Connector.hpp"
#include "Rest.hpp"
#include "CheckNoOverlaps.hpp"
#include "OutputPipe.hpp"

using namespace bloomrepeats;

class RestPipe : public Pipe {
public:
    RestPipe() {
        add(new AddBlocks, "target=other");
        add(new Connector, "target=other");
        add(new Rest);
        add(new CheckNoOverlaps);
        add(new OutputPipe, "--skip-rest:=1");
    }
};

int main(int argc, char** argv) {
    return process(argc, argv, new RestPipe, "Select unique fragments");
}

