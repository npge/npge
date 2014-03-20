/*
 * bloomrepeats, Find genomic repeats, using Bloom filter based prefiltration
 * Copyright (C) 2013 Boris Nagaev
 *
 * See the LICENSE file for terms of use.
 */

#include <boost/filesystem.hpp>

#include "FileCopy.hpp"
#include "name_to_stream.hpp"

namespace bloomrepeats {

FileCopy::FileCopy() {
    add_opt("src", "Source file to copy", std::string(), true);
    add_opt("dst", "Destination of copy", std::string(":null"));
}

bool FileCopy::run_impl() const {
    std::string src = opt_value("src").as<std::string>();
    std::string dst = opt_value("dst").as<std::string>();
    if (dst != ":null") {
        using namespace boost::filesystem;
        copy_file(src, dst);
    }
    return false;
}

}

