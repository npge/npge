/*
 * bloomrepeats, Find genomic repeats, using Bloom filter based prefiltration
 * Copyright (C) 2014 Boris Nagaev
 *
 * See the LICENSE file for terms of use.
 */

#include <boost/foreach.hpp>
#include <boost/scoped_ptr.hpp>

#include "PrintBlockSetAlignment.hpp"
#include "block_set_alignment.hpp"
#include "BlockSet.hpp"

namespace bloomrepeats {

PrintBlockSetAlignment::PrintBlockSetAlignment():
    file_writer_(this, "out-bs-aln",
                 "Output file with block set alignment") {
    add_opt("bs-aln-chr", "Chromosome used for block set alignment",
            std::string("chr1"));
    add_opt("bs-aln-blocks", "Print block names in alignment "
            "(else fragments)", true);
}

bool PrintBlockSetAlignment::run_impl() const {
    std::ostream& out = file_writer_.output();
    std::string chr = opt_value("bs-aln-chr").as<std::string>();
    bool blocks = opt_value("bs-aln-blocks").as<bool>();
    bsa_print(out, block_set()->bsa(chr), blocks);
    return false;
}

const char* PrintBlockSetAlignment::name_impl() const {
    return "Print block set alignment";
}

}

