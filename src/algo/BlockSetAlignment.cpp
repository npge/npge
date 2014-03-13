/*
 * bloomrepeats, Find genomic repeats, using Bloom filter based prefiltration
 * Copyright (C) 2014 Boris Nagaev
 *
 * See the LICENSE file for terms of use.
 */

#include <boost/foreach.hpp>

#include "BlockSetAlignment.hpp"
#include "block_set_alignment.hpp"
#include "Fragment.hpp"
#include "Block.hpp"
#include "Sequence.hpp"

namespace bloomrepeats {

BlockSetAlignment::BlockSetAlignment():
    file_writer_(this, "out-bs-aln",
                 "Output file with block set alignment") {
    add_opt("bs-aln-chr", "Chromosome used for block set alignment",
            std::string("chr1"));
}

bool BlockSetAlignment::run_impl() const {
    std::string chr = opt_value("bs-aln-chr").as<std::string>();
    std::ostream& out = file_writer_.output();
    BSA rows, aln;
    bsa_make_rows(rows, *block_set(), chr);
    bsa_make_aln(aln, rows);
    BOOST_FOREACH (const BSA::value_type& seq_and_row, aln) {
        Sequence* seq = seq_and_row.first;
        const BSRow& row = seq_and_row.second;
        out << ((row.ori == 1) ? '+' : '-');
        out << seq->genome();
        BOOST_FOREACH (Fragment* fragment, row.fragments) {
            out << '\t';
            if (fragment) {
                out << fragment->block()->name();
            } else {
                out << '-';
            }
        }
        out << "\n";
    }
    return false;
}

}

