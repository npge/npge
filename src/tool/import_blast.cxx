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
#include "ImportBlastHits.hpp"
#include "OutputPipe.hpp"

using namespace bloomrepeats;

class ImportBlastHitsPipe : public Pipe {
public:
    ImportBlastHitsPipe() {
        add(new AddSequences, "target=other");
        AddBlocks* ab = new AddBlocks(/* keep_alignment */ true);
        ab->add_ignored_option("import-alignment");
        add(ab, "target=other");
        add(new ImportBlastHits);
        add(new OutputPipe);
    }
};

int main(int argc, char** argv) {
    return process(argc, argv, new ImportBlastHitsPipe,
                   "Import blast hits");
}

