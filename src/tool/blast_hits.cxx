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
#include "AddBlastBlocks.hpp"
#include "UniqueNames.hpp"
#include "Output.hpp"

using namespace bloomrepeats;

class BlastHitsPipe : public Pipe {
public:
    BlastHitsPipe() {
        set_empty_block_set();
        set_empty_other();
        add(new AddSequences, OTHER_TO_THIS);
        add(new AddBlocks, OTHER_TO_THIS);
        add(new AddBlastBlocks);
        add(new UniqueNames);
        add(new Output);
    }
};

int main(int argc, char** argv) {
    return process(argc, argv, new BlastHitsPipe,
                   "Print blocks found by blast");
}

