/*
 * NPG-explorer, Nucleotide PanGenome explorer
 * Copyright (C) 2014 Boris Nagaev
 *
 * See the LICENSE file for terms of use.
 */

#include <cctype>
#include <streambuf>
#include <algorithm>
#include <boost/algorithm/string/trim.hpp>

#include "FastaReader.hpp"

namespace npge {

FastaReader::FastaReader(std::istream& input):
    input_(input) {
}

bool FastaReader::read_one_sequence() {
    bool in_sequence = false;
    std::streampos line_start = input_.tellg();
    for (std::string line; std::getline(input_, line);) {
        boost::algorithm::trim(line);
        std::streamoff line_size = line.size();
        if (line.empty()) {
            empty_line_found(); // TODO test for this
            found_empty_line_ = true;
        } else if (line_size >= 1 && line[0] == '>') {
            if (!in_sequence) {
                in_sequence = true;
                std::string name;
                std::string description;
                if (line.size() >= 2) {
                    size_t sp = std::string::npos;
                    size_t s = line.size();
                    for (size_t i = 0; i < s; i++) {
                        if (isspace(line[i])) {
                            sp = i;
                            break;
                        }
                    }
                    name = line.substr(1, sp - 1);
                    if (sp != std::string::npos) {
                        size_t dp = sp + 1;
                        while (dp < s && isspace(line[dp])) {
                            dp += 1;
                        }
                        if (dp < s) {
                            description = line.substr(dp);
                        }
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

bool FastaReader::read_until_empty_line() {
    bool result = false;
    while (true) {
        found_empty_line_ = false;
        bool ok = read_one_sequence();
        result |= ok;
        if (found_empty_line_ || !ok) {
            break;
        }
    }
    return result;
}

bool FastaReader::read_all_sequences() {
    bool result = false;
    while (read_one_sequence()) {
        result = true;
    }
    return result;
}

void FastaReader::empty_line_found() {
}

}

