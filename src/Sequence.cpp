/*
 * bloomrepeats, Find genomic repeats, using Bloom filter based prefiltration
 * Copyright (C) 2012 Boris Nagaev
 *
 * See the LICENSE file for terms of use.
 */

#include <fstream>
#include <streambuf>
#include <boost/algorithm/string/replace.hpp>
#include <boost/algorithm/string/trim.hpp>

#include "Sequence.hpp"

namespace bloomrepeats {

InMemorySequence::InMemorySequence(const std::string& filename) {
    std::ifstream file(filename.c_str());
    bool in_sequence = false;
    for (std::string line; std::getline(file, line);) {
        boost::algorithm::trim(line);
        boost::algorithm::replace_all(line, " ", "");
        boost::algorithm::replace_all(line, "\t", "");
        if (line[0] == '>') {
            if (data_.empty()) {
                in_sequence = true;
            } else {
                break;
            }
        } else if (in_sequence) {
            data_ += line;
        }
    }
}

const char* InMemorySequence::get(size_t start, size_t& length) const {
    if (data_.size() < start + length) {
        length = data_.size() - start;
    }
    return data_.c_str() + start;
}

size_t InMemorySequence::approximate_size() const {
    return data_.size();
}

}

