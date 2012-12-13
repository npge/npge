/*
 * bloomrepeats, Find genomic repeats, using Bloom filter based prefiltration
 * Copyright (C) 2012 Boris Nagaev
 *
 * See the LICENSE file for terms of use.
 */

#include <cstdio>

#include "FileWriter.hpp"
#include "temp_file.hpp"

namespace bloomrepeats {

FileWriter::FileWriter() {
    set_rand_name();
    set_remove_after(true);
}

FileWriter::~FileWriter() {
    if (get_remove_after()) {
        remove_file();
    }
}

void FileWriter::set_output_file(const std::string& output_file,
                                 bool remove_prev) {
    if (remove_prev) {
        remove_file();
    }
    output_file_ = output_file;
}

void FileWriter::set_rand_name(bool remove_prev) {
    set_output_file(temp_file(), remove_prev);
}

void FileWriter::remove_file() {
    remove(output_file().c_str());
}

void FileWriter::set_remove_after(bool value) {
    remove_after_ = value;
}

}

