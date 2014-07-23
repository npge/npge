/*
 * NPG-explorer, Nucleotide PanGenome explorer
 * Copyright (C) 2013 Boris Nagaev
 *
 * See the LICENSE file for terms of use.
 */

#include "FileCopy.hpp"
#include "name_to_stream.hpp"

namespace npge {

FileCopy::FileCopy() {
    add_opt("src", "Source file to copy", std::string(), true);
    add_opt("dst", "Destination of copy", std::string(":null"));
}

void FileCopy::run_impl() const {
    std::string src = opt_value("src").as<std::string>();
    std::string dst = opt_value("dst").as<std::string>();
    if (dst != ":null") {
        copy_file(src, dst);
    }
}

const char* FileCopy::name_impl() const {
    return "Copy file";
}

}

