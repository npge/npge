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
#include <boost/algorithm/string/case_conv.hpp>

#include "Sequence.hpp"
#include "Fragment.hpp"

namespace bloomrepeats {

void Sequence::make_first_fragment(Fragment& f, size_t fragment_size,
                                   int only_ori) const {
    f.set_min_pos(0 - 1);
    f.set_max_pos(fragment_size - 1 - 1);
    f.set_ori(only_ori ? : 1);
}

bool Sequence::next_fragment(Fragment& f) const {
    f.set_ori(-f.ori());
    if (f.ori() == -1) {
        f.set_min_pos(f.min_pos() + 1);
        f.set_max_pos(f.max_pos() + 1);
    }
    return f.max_pos() < size();
}

bool Sequence::next_fragment_keeping_ori(Fragment& f) const {
    f.set_min_pos(f.min_pos() + 1);
    f.set_max_pos(f.max_pos() + 1);
    return f.max_pos() < size();
}

InMemorySequence::InMemorySequence(const std::string& filename, int) {
    std::ifstream file(filename.c_str());
    read_from_file(file);
}

InMemorySequence::InMemorySequence(std::istream& input) {
    read_from_file(input);
}

InMemorySequence::InMemorySequence(const std::string& data):
    data_(boost::algorithm::to_lower_copy(data)) {
    set_size(data_.size());
}

const char* InMemorySequence::get(size_t start, size_t& length) const {
    if (data_.size() < start + length) {
        length = data_.size() - start;
    }
    return data_.c_str() + start;
}

void InMemorySequence::read_from_file(std::istream& input) {
    bool in_sequence = false;
    for (std::string line; std::getline(input, line);) {
        std::streamoff line_size = line.size();
        boost::algorithm::trim(line);
        boost::algorithm::replace_all(line, " ", "");
        boost::algorithm::replace_all(line, "\t", "");
        if (line[0] == '>') {
            if (data_.empty()) {
                in_sequence = true;
            } else {
                // go to the beginning of current line
                input.seekg(input.tellg() - line_size - std::streamoff(2));
                break;
            }
        } else if (in_sequence) {
            boost::algorithm::to_lower(line);
            data_ += line;
        }
    }
    set_size(data_.size());
}

}

