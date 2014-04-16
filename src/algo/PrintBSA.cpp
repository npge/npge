/*
 * bloomrepeats, Find genomic repeats, using Bloom filter based prefiltration
 * Copyright (C) 2014 Boris Nagaev
 *
 * See the LICENSE file for terms of use.
 */

#include <boost/foreach.hpp>
#include <boost/scoped_ptr.hpp>

#include "PrintBSA.hpp"
#include "bsa_algo.hpp"
#include "BlockSet.hpp"

namespace bloomrepeats {

PrintBSA::PrintBSA():
    file_writer_(this, "out-bsa",
                 "Output file with block set alignment") {
    add_opt("bsa-blocks", "Print block names in alignment "
            "(else fragments)", true);
    declare_bs("target", "Target blockset");
}

void PrintBSA::run_impl() const {
    std::ostream& out = file_writer_.output();
    bool blocks = opt_value("bsa-blocks").as<bool>();
    BOOST_FOREACH (std::string bsa_name, block_set()->bsas()) {
        bsa_print(out, block_set()->bsa(bsa_name), bsa_name, blocks);
    }
}

const char* PrintBSA::name_impl() const {
    return "Print block set alignment";
}

}

