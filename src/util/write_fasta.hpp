/*
 * bloomrepeats, Find genomic repeats, using Bloom filter based prefiltration
 * Copyright (C) 2014 Boris Nagaev
 *
 * See the LICENSE file for terms of use.
 */

#ifndef BR_WRITE_FASTA_HPP_
#define BR_WRITE_FASTA_HPP_

#include <iosfwd>
#include <string>

namespace bloomrepeats {

void write_fasta(std::ostream& out, const std::string& name,
                 const std::string& description,
                 const std::string& text,
                 int line = 0);

}

#endif

