/*
 * bloomrepeats, Find genomic repeats, using Bloom filter based prefiltration
 * Copyright (C) 2012 Boris Nagaev
 *
 * See the LICENSE file for terms of use.
 */

#include <fstream>
#include <iostream>
#include <streambuf>

#include "read_file.hpp"

namespace bloomrepeats {

std::string read_stream(std::istream& stream) {
    return std::string((std::istreambuf_iterator<char>(stream)),
                       std::istreambuf_iterator<char>());
}

std::string read_file(const std::string& filename) {
    std::ifstream stream(filename.c_str());
    return read_stream(stream);
}

std::string read_stdin() {
    return read_stream(std::cin);
}

}

