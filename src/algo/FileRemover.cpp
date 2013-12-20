/*
 * bloomrepeats, Find genomic repeats, using Bloom filter based prefiltration
 * Copyright (C) 2013 Boris Nagaev
 *
 * See the LICENSE file for terms of use.
 */

#include "FileRemover.hpp"
#include "name_to_stream.hpp"

namespace bloomrepeats {

FileRemover::FileRemover() {
    add_opt("filename", "File to be removed", std::string(), true);
}

bool FileRemover::run_impl() const {
    remove_file(opt_value("filename").as<std::string>());
    return false;
}

}

