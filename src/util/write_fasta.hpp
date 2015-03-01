/*
 * NPG-explorer, Nucleotide PanGenome explorer
 * Copyright (C) 2012-2015 Boris Nagaev
 *
 * See the LICENSE file for terms of use.
 */

#ifndef NPGE_WRITE_FASTA_HPP_
#define NPGE_WRITE_FASTA_HPP_

#include <iosfwd>
#include <string>

namespace npge {

void write_fasta(std::ostream& out, const std::string& name,
                 const std::string& description,
                 const std::string& text,
                 int line = 0);

}

#endif

