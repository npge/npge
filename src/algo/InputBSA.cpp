/*
 * NPG-explorer, Nucleotide PanGenome explorer
 * Copyright (C) 2012-2015 Boris Nagaev
 *
 * See the LICENSE file for terms of use.
 */

#include <boost/foreach.hpp>

#include "InputBSA.hpp"
#include "bsa_algo.hpp"
#include "BlockSet.hpp"

namespace npge {

InputBSA::InputBSA():
    file_reader_(this, "in-bsa", "input file(s) with block "
                 "set alignments") {
    declare_bs("target", "Target blockset");
}

void InputBSA::run_impl() const {
    BOOST_FOREACH (std::istream& input_file, file_reader_) {
        bsa_input(*block_set(), input_file);
    }
}

const char* InputBSA::name_impl() const {
    return "Input blockset alignment";
}

}

