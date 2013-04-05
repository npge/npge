/*
 * bloomrepeats, Find genomic repeats, using Bloom filter based prefiltration
 * Copyright (C) 2012 Boris Nagaev
 *
 * See the LICENSE file for terms of use.
 */

#include "process.hpp"
#include "Pipe.hpp"
#include "AddBlocks.hpp"
#include "Filter.hpp"
#include "OutputPipe.hpp"

using namespace bloomrepeats;

class FilterPipe : public Pipe {
public:
    FilterPipe() {
        add(new AddBlocks);
        add(new Filter);
        add(new OutputPipe);
    }
};

int main(int argc, char** argv) {
    return process(argc, argv, new FilterPipe,
                   "Filter blocks");
}

