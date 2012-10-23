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
#include "Swap.hpp"
#include "ImportBlastHits.hpp"
#include "UniqueNames.hpp"
#include "Output.hpp"

using namespace bloomrepeats;

class FromBlastPipe : public Pipe {
public:
    FromBlastPipe() {
        BlockSetPtr reference = boost::make_shared<BlockSet>();
        set_block_set(reference);
        add(new AddSequences);
        add(new AddBlocks(/* keep_alignment */ true));
        BlockSetPtr blast_hits = boost::make_shared<BlockSet>();
        add(new Swap(blast_hits));
        blast_hits.swap(reference);
        add(new ImportBlastHits(reference));
        add(new UniqueNames);
        add(new Output);
    }
};

int main(int argc, char** argv) {
    process(argc, argv, new FromBlastPipe, "Import blast hits");
}

