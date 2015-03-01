/*
 * NPG-explorer, Nucleotide PanGenome explorer
 * Copyright (C) 2012-2015 Boris Nagaev
 *
 * See the LICENSE file for terms of use.
 */

#ifndef NPGE_IMPORT_BLAST_HITS_HPP_
#define NPGE_IMPORT_BLAST_HITS_HPP_

#include "global.hpp"
#include "Processor.hpp"
#include "FileReader.hpp"

namespace npge {

/** Add blocks from blast hits (blast output -m 8).
\note This processor depends on processor Read.
*/
class ImportBlastHits : public Processor {
public:
    /** Constructor */
    ImportBlastHits();

protected:
    void run_impl() const;

    const char* name_impl() const;

private:
    FileReader file_reader_;
};

}

#endif

