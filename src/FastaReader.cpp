/*
 * bloomrepeats, Find genomic repeats, using Bloom filter based prefiltration
 * Copyright (C) 2012 Boris Nagaev
 *
 * See the LICENSE file for terms of use.
 */

#include <cctype>
#include <fstream>
#include <streambuf>
#include <algorithm>
#include <boost/algorithm/string/trim.hpp>

#include "FastaReader.hpp"

namespace bloomrepeats {

FastaReader::FastaReader(std::istream& input):
    input_(input)
{ }

bool FastaReader::read_one_sequence() {
    bool in_sequence = false;
    std::streampos line_start = input_.tellg();
    for (std::string line; std::getline(input_, line);) {
        boost::algorithm::trim(line);
        std::streamoff line_size = line.size();
        if (line.empty()) {
            empty_line_found(); // TODO test for this
        } else if (line_size >= 1 && line[0] == '>') {
            if (!in_sequence) {
                in_sequence = true;
                std::string name;
                std::string description;
                if (line.size() >= 2) {
                    size_t sp = line.find(' ');
                    name = line.substr(1, sp - 1);
                    if (sp != std::string::npos && sp + 1 < line.size()) {
                        description = line.substr(sp + 1);
                    }
                }
                new_sequence(name, description);
            } else {
                // go to the beginning of current line
                input_.seekg(line_start);
                return true;
            }
        } else if (in_sequence) {
            line.erase(std::remove_if(line.begin(), line.end(), isspace),
                       line.end());
            grow_sequence(line);
        }
        line_start = input_.tellg();
    }
    return in_sequence;
}

bool FastaReader::read_all_sequences() {
    bool result = false;
    while (read_one_sequence()) {
        result = true;
    }
    return result;
}

void FastaReader::empty_line_found()
{ }

}

