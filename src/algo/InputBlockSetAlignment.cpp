/*
 * bloomrepeats, Find genomic repeats, using Bloom filter based prefiltration
 * Copyright (C) 2014 Boris Nagaev
 *
 * See the LICENSE file for terms of use.
 */

#include <boost/foreach.hpp>

#include "InputBlockSetAlignment.hpp"
#include "block_set_alignment.hpp"
#include "BlockSet.hpp"

namespace bloomrepeats {

InputBlockSetAlignment::InputBlockSetAlignment():
    file_reader_(this, "in-bsa", "input file(s) with block "
                 "set alignments")
{ }

bool InputBlockSetAlignment::run_impl() const {
    BOOST_FOREACH (std::istream& input_file, file_reader_) {
        bsa_input(*block_set(), input_file);
    }
    return false;
}

const char* InputBlockSetAlignment::name_impl() const {
    return "Input block set alignment";
}

}

