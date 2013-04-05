/*
 * bloomrepeats, Find genomic repeats, using Bloom filter based prefiltration
 * Copyright (C) 2012 Boris Nagaev
 *
 * See the LICENSE file for terms of use.
 */

#include "process.hpp"
#include "Pipe.hpp"
#include "AddBlocks.hpp"
#include "ExternalAligner.hpp"
#include "Output.hpp"

using namespace bloomrepeats;

class ExternalAlignerPipe : public Pipe {
public:
    ExternalAlignerPipe() {
        add(new AddBlocks);
        add(new ExternalAligner);
        add(new Output, "--export-alignment:=1");
    }
};

int main(int argc, char** argv) {
    return process(argc, argv, new ExternalAlignerPipe,
                   "Align all blocks (by default, with fast mafft)");
}

