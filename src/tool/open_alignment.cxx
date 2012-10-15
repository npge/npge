/*
 * bloomrepeats, Find genomic repeats, using Bloom filter based prefiltration
 * Copyright (C) 2012 Boris Nagaev
 *
 * See the LICENSE file for terms of use.
 */

#include "process.hpp"
#include "BlockSet.hpp"
#include "Alignment.hpp"
#include "Pipe.hpp"
#include "AddSequences.hpp"
#include "AddBlocks.hpp"
#include "Output.hpp"

using namespace bloomrepeats;

class OpenAlignmentPipe : public Pipe {
public:
    OpenAlignmentPipe() {
        set_block_set(boost::make_shared<BlockSet>());
        add(new AddSequences);
        AlignmentPtr alignment = boost::make_shared<Alignment>();
        add(new AddBlocks(alignment));
        add(new Output);
    }
};

int main(int argc, char** argv) {
    return process(argc, argv, new OpenAlignmentPipe,
                   "Read alignment and write block set (test)");
}

