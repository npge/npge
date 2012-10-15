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
#include "Connector.hpp"
#include "Rest.hpp"
#include "Union.hpp"
#include "Swap.hpp"
#include "UniqueNames.hpp"
#include "CheckNoOverlaps.hpp"
#include "Output.hpp"

using namespace bloomrepeats;

class RestPipe : public Pipe {
public:
    RestPipe() {
        set_block_set(boost::make_shared<BlockSet>());
        add(new AddSequences);
        add(new AddBlocks);
        add(new Connector);
        BlockSetPtr pre_rest = boost::make_shared<BlockSet>();
        add(new Swap(pre_rest));
        add(new Rest(pre_rest));
        add(new Union(pre_rest));
        add(new CheckNoOverlaps);
        add(new UniqueNames);
        add(new Output);
    }
};

int main(int argc, char** argv) {
    return process(argc, argv, new RestPipe,
                   "Add unique fragments to current block set");
}

