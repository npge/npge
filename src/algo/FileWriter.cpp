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

void FileWriter::set_file(const std::string& file, bool remove_prev) {
    if (remove_prev) {
        remove_file();
    }
    file_ = file;
}

void FileWriter::set_rand_name(bool remove_prev) {
    set_file(temp_file(), remove_prev);
}

void FileWriter::remove_file() {
    remove(file().c_str());
}

void FileWriter::set_remove_after(bool value) {
    remove_after_ = value;
}

}

