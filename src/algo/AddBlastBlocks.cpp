/*
 * bloomrepeats, Find genomic repeats, using Bloom filter based prefiltration
 * Copyright (C) 2012 Boris Nagaev
 *
 * See the LICENSE file for terms of use.
 */

#include "AddBlastBlocks.hpp"
#include "Consensus.hpp"
#include "BlastRunner.hpp"
#include "ImportBlastHits.hpp"
#include "BlockSet.hpp"

namespace bloomrepeats {

AddBlastBlocks::AddBlastBlocks() {
    Consensus* consensus = new Consensus;
    consensus->set_no_options(true);
    add(consensus);
    BlastRunner* blast_runner = new BlastRunner;
    blast_runner->set_input_file(consensus->file());
    blast_runner->set_no_options(true);
    add(blast_runner);
    ImportBlastHits* import_blast = new ImportBlastHits;
    import_blast->set_processor(this);
    import_blast->set_input_file(blast_runner->file());
    import_blast->set_no_options(true);
    add(import_blast);
}

}

