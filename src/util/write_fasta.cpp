/*
 * NPG-explorer, Nucleotide PanGenome explorer
 * Copyright (C) 2012-2015 Boris Nagaev
 *
 * See the LICENSE file for terms of use.
 */

#include <ostream>

namespace npge {

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
            size_t l = (i == lines - 1) ? std::string::npos : line;
            out << text.substr(i * line, l) << '\n';
        }
    }
}

}

