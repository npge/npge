/*
 * NPG-explorer, Nucleotide PanGenome explorer
 * Copyright (C) 2012 Boris Nagaev
 *
 * See the LICENSE file for terms of use.
 */

#ifndef NPGE_ADD_GENES_HPP_
#define NPGE_ADD_GENES_HPP_

#include "Processor.hpp"
#include "FileReader.hpp"

namespace npge {

/** Add genes from EBI genes description.
Sequence accession numbers are taken from Sequence.ac().
*/
class AddGenes : public Processor {
public:
    /** Constructor */
    AddGenes();

protected:
    void run_impl() const;

    const char* name_impl() const;

private:
    FileReader file_reader_;
};

}

#endif

