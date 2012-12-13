/*
 * bloomrepeats, Find genomic repeats, using Bloom filter based prefiltration
 * Copyright (C) 2012 Boris Nagaev
 *
 * See the LICENSE file for terms of use.
 */

#include "FileReader.hpp"

namespace bloomrepeats {

void FileReader::set_input_file(const std::string& input_file) {
    input_files_.clear();
    input_files_.push_back(input_file);
}

}

