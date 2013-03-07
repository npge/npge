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
#include "BlockInfo.hpp"

using namespace bloomrepeats;

class BlockInfoPipe : public Pipe {
public:
    BlockInfoPipe() {
        set_empty_block_set();
        add(new AddSequences);
        add(new AddBlocks);
        add(new BlockInfo);
    }
};

int main(int argc, char** argv) {
    return process(argc, argv, new BlockInfoPipe,
                   "Print information about blocks");
}

