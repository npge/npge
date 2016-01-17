/*
 * NPG-explorer, Nucleotide PanGenome explorer
 * Copyright (C) 2012-2016 Boris Nagaev
 *
 * See the LICENSE file for terms of use.
 */

#ifndef NPGE_FASTA_READER_HPP_
#define NPGE_FASTA_READER_HPP_

#include <iosfwd>
#include <string>

#include "global.hpp"

namespace npge {

class FastaReader {
public:
    FastaReader(std::istream& input);

    bool read_one_sequence();

    bool read_until_empty_line();

    bool read_all_sequences();

protected:
    virtual void new_sequence(const std::string& name,
                              const std::string& description) = 0;

    virtual void grow_sequence(const std::string& data) = 0;

    virtual void empty_line_found();

private:
    std::istream& input_;
    bool found_empty_line_;
};

}

#endif

