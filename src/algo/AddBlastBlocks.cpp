/*
 * bloomrepeats, Find genomic repeats, using Bloom filter based prefiltration
 * Copyright (C) 2012 Boris Nagaev
 *
 * See the LICENSE file for terms of use.
 */

#include "SequencesFromOther.hpp"
#include "AddBlastBlocks.hpp"
#include "Consensus.hpp"
#include "BlastRunner.hpp"
#include "ImportBlastHits.hpp"
#include "UniqueNames.hpp"
#include "BlockSet.hpp"

namespace bloomrepeats {

AddBlastBlocks::AddBlastBlocks(BlockSetPtr source):
        Pipe(source) {
    add(new UniqueNames, OTHER_TO_THIS);
    Consensus* consensus = new Consensus;
    consensus->set_no_options(true);
    add(consensus, OTHER_TO_THIS);
    BlastRunner* blast_runner = new BlastRunner;
    blast_runner->set_input_file(consensus->output_file());
    blast_runner->set_no_options(true);
    add(blast_runner, OTHER_TO_THIS);
    add(new SequencesFromOther(BlockSetPtr()));
    ImportBlastHits* import_blast = new ImportBlastHits;
    import_blast->set_input_file(blast_runner->output_file());
    import_blast->set_no_options(true);
    add(import_blast);
}

}

