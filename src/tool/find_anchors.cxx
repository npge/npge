/*
 * bloomrepeats, Find genomic repeats, using Bloom filter based prefiltration
 * Copyright (C) 2012 Boris Nagaev
 *
 * See the LICENSE file for terms of use.
 */

#include "process.hpp"
#include "Pipe.hpp"
#include "AddSequences.hpp"
#include "AnchorFinder.hpp"
#include "OutputPipe.hpp"

using namespace bloomrepeats;

class AnchorsPipe : public Pipe {
public:
    AnchorsPipe() {
        set_empty_block_set();
        add(new AddSequences);
        add(new AnchorFinder);
        add(new OutputPipe);
    }
};

int main(int argc, char** argv) {
    return process(argc, argv, new AnchorsPipe, "Find anchors", "in-seqs");
}

