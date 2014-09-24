/*
 * NPG-explorer, Nucleotide PanGenome explorer
 * Copyright (C) 2014 Boris Nagaev
 *
 * See the LICENSE file for terms of use.
 */

#include <boost/foreach.hpp>
#include <boost/scoped_ptr.hpp>

#include "PrintBSA.hpp"
#include "bsa_algo.hpp"
#include "BlockSet.hpp"

namespace npge {

PrintBSA::PrintBSA():
    file_writer_(this, "out-bsa",
                 "Output file with blockset alignment") {
    add_opt("bsa-blocks", "Print block names in alignment "
            "(else fragments)", true);
    add_opt("bsa-conservative", "Print conservative columns line",
            true);
    add_opt("bsa-orientation", "Print orientation after fragment",
            true);
    declare_bs("target", "Target blockset");
}

void PrintBSA::run_impl() const {
    std::ostream& out = file_writer_.output();
    bool blocks = opt_value("bsa-blocks").as<bool>();
    bool conservative = opt_value("bsa-conservative").as<bool>();
    bool orientation = opt_value("bsa-conservative").as<bool>();
    BOOST_FOREACH (std::string bsa_name, block_set()->bsas()) {
        const BSA& bsa = block_set()->bsa(bsa_name);
        bsa_print(out, bsa, bsa_name, blocks, orientation);
        if (conservative) {
            bsa_print_conservative(out, bsa, bsa_name);
        }
    }
}

const char* PrintBSA::name_impl() const {
    return "Print blockset alignment";
}

}

