/*
 * NPG-explorer, Nucleotide PanGenome explorer
 * Copyright (C) 2013 Boris Nagaev
 *
 * See the LICENSE file for terms of use.
 */

#include "FileRemover.hpp"
#include "name_to_stream.hpp"

namespace npge {

FileRemover::FileRemover() {
    add_opt("filename", "File to be removed", std::string(), true);
}

void FileRemover::run_impl() const {
    remove_file(opt_value("filename").as<std::string>());
}

const char* FileRemover::name_impl() const {
    return "Remove file";
}

}

