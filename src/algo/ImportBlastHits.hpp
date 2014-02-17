/*
 * bloomrepeats, Find genomic repeats, using Bloom filter based prefiltration
 * Copyright (C) 2012 Boris Nagaev
 *
 * See the LICENSE file for terms of use.
 */

#ifndef BR_IMPORT_BLAST_HITS_HPP_
#define BR_IMPORT_BLAST_HITS_HPP_

#include "global.hpp"
#include "Processor.hpp"
#include "FileReader.hpp"

namespace bloomrepeats {

/** Add blocks from blast hits (blast output -m 8).
\note This processor depends on AddBlocks(keep_alignment = true).
*/
class ImportBlastHits : public Processor {
public:
    /** Constructor.
    \param block_set The block set, passed to blast.
    \param min_length Min accepted length of blast hit
    \param min_ident Min accepted identity of blast hit
    \param max_evalue Max accepted e-value of blast hit
    */
    ImportBlastHits(const BlockSetPtr& block_set = BlockSetPtr(),
                    int min_length = 100,
                    double min_ident = 0.9, double max_evalue = 0.001);

protected:
    bool run_impl() const;

    const char* name_impl() const;

private:
    FileReader file_reader_;
};

}

#endif

