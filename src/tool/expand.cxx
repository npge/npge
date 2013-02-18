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
#include "Connector.hpp"
#include "FragmentsExpander.hpp"
#include "UniqueNames.hpp"
#include "Output.hpp"

using namespace bloomrepeats;

class FragmentsExpanderPipe : public Pipe {
public:
    FragmentsExpanderPipe() {
        set_empty_block_set();
        add(new AddSequences);
        add(new AddBlocks);
        add(new Connector);
        add(new FragmentsExpander);
        add(new UniqueNames);
        add(new Output);
    }
};

int main(int argc, char** argv) {
    return process(argc, argv, new FragmentsExpanderPipe, "Expand fragments");
}

