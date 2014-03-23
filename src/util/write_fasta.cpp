/*
 * bloomrepeats, Find genomic repeats, using Bloom filter based prefiltration
 * Copyright (C) 2014 Boris Nagaev
 *
 * See the LICENSE file for terms of use.
 */

#include <ostream>

namespace bloomrepeats {

void write_fasta(std::ostream& out, const std::string& name,
                 const std::string& description,
                 const std::string& text,
                 int line) {
    out << '>' << name << ' ' << description << '\n';
    if (line == 0) {
        out << text << '\n';
    } else {
        int lines = (text.length() + line - 1) / line;
        for (int i = 0; i < lines; i++) {
            int l = (i == lines - 1) ? std::string::npos : line;
            out << text.substr(i * line, l) << '\n';
        }
    }
}

}

