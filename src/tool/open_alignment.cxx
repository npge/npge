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
#include "OutputPipe.hpp"

using namespace bloomrepeats;

class OpenAlignmentPipe : public Pipe {
public:
    OpenAlignmentPipe() {
        set_empty_block_set();
        add(new AddSequences);
        add(new AddBlocks(/* keep_alignment */ true));
        add(new OutputPipe);
    }
};

int main(int argc, char** argv) {
    process(argc, argv, new OpenAlignmentPipe,
            "Read alignment and write block set (test)");
}

