/*
 * bloomrepeats, Find genomic repeats, using Bloom filter based prefiltration
 * Copyright (C) 2012 Boris Nagaev
 *
 * See the LICENSE file for terms of use.
 */

#include "process.hpp"
#include "Pipe.hpp"
#include "AddBlocks.hpp"
#include "AddBlastBlocks.hpp"
#include "CleanUp.hpp"
#include "CheckNoOverlaps.hpp"
#include "OutputPipe.hpp"

using namespace bloomrepeats;

class AddBlastPipe : public Pipe {
public:
    AddBlastPipe() {
        add(new AddBlocks, "--import-alignment:=1");
        add(new AddBlastBlocks, "other=target target=target");
        add(new CleanUp);
        add(new CheckNoOverlaps);
        add(new OutputPipe);
    }
};

int main(int argc, char** argv) {
    return process(argc, argv, new AddBlastPipe,
                   "Add blocks found by blast");
}

