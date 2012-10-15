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
#include "AnchorFinder.hpp"
#include "CleanUp.hpp"
#include "CheckNoOverlaps.hpp"
#include "UniqueNames.hpp"
#include "Output.hpp"

using namespace bloomrepeats;

class BloomRepeatsPipe : public Pipe {
public:
    BloomRepeatsPipe() {
        set_block_set(boost::make_shared<BlockSet>());
        add(new AddSequences);
        add(new AnchorFinder);
        add(new CleanUp);
        add(new CheckNoOverlaps);
        add(new UniqueNames);
        add(new Output);
    }
};

int main(int argc, char** argv) {
    return process(argc, argv, new BloomRepeatsPipe,
                   "Find and expand anchors", "in-seqs");
}

