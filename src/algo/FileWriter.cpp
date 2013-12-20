/*
 * bloomrepeats, Find genomic repeats, using Bloom filter based prefiltration
 * Copyright (C) 2012 Boris Nagaev
 *
 * See the LICENSE file for terms of use.
 */

#include <boost/filesystem.hpp>

#include "FileWriter.hpp"
#include "name_to_stream.hpp"
#include "temp_file.hpp"

namespace bloomrepeats {

FileWriter::FileWriter() {
    set_remove_after(true);
}

FileWriter::~FileWriter() {
    if (get_remove_after()) {
        remove_file(output_file());
    }
}

void FileWriter::set_output_file(const std::string& output_file,
                                 bool remove_prev) {
    if (get_remove_after()) {
        remove_file(this->output_file());
    }
    output_file_ = output_file;
    output_.reset();
}

void FileWriter::set_rand_name(bool remove_prev) {
    set_output_file(temp_file(), remove_prev);
}

void FileWriter::set_remove_after(bool value) {
    remove_after_ = value;
}

std::ostream& FileWriter::output() const {
    using namespace boost::filesystem;
    if (!output_ || (output_file()[0] != ':' && !exists(output_file()))) {
        output_ = name_to_ostream(output_file());
    }
    return *output_;
}

}

