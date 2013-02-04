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
#include "AddBlastBlocks.hpp"
#include "CleanUp.hpp"
#include "UniqueNames.hpp"
#include "CheckNoOverlaps.hpp"
#include "Output.hpp"

using namespace bloomrepeats;

class AddBlastPipe : public Pipe {
public:
    AddBlastPipe() {
        set_empty_block_set();
        set_empty_other();
        add(new AddSequences);
        add(new AddBlocks);
        add(new AddBlastBlocks, THIS_TO_OTHER | THIS_TO_THIS);
        add(new CleanUp);
        add(new CheckNoOverlaps);
        add(new UniqueNames);
        add(new Output);
    }
};

int main(int argc, char** argv) {
    return process(argc, argv, new AddBlastPipe,
                   "Add blocks found by blast");
}

