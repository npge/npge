/*
 * NPG-explorer, Nucleotide PanGenome explorer
 * Copyright (C) 2014 Boris Nagaev
 *
 * See the LICENSE file for terms of use.
 */

#ifndef BR_WRITE_FASTA_HPP_
#define BR_WRITE_FASTA_HPP_

#include <iosfwd>
#include <string>

namespace npge {

void write_fasta(std::ostream& out, const std::string& name,
                 const std::string& description,
                 const std::string& text,
                 int line = 0);

}

#endif

