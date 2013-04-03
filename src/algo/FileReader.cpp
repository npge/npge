/*
 * bloomrepeats, Find genomic repeats, using Bloom filter based prefiltration
 * Copyright (C) 2012 Boris Nagaev
 *
 * See the LICENSE file for terms of use.
 */

#include <fstream>

#include "FileReader.hpp"
#include "Exception.hpp"

namespace bloomrepeats {

typedef FileReader::const_iterator FRCI;

FRCI& FRCI::operator++() {
    stream_.reset();
    index_++;
    return *this;
}

FRCI& FRCI::operator++(int) {
    return ++(*this);
}

bool FRCI::operator==(const FRCI& other) {
    return reader_ == other.reader_ && index_ == other.index_;
}

bool FRCI::operator!=(const FRCI& other) {
    return !(*this == other);
}

std::istream& FRCI::operator*() {
    if (!stream_) {
        stream_.reset(reader_->input_stream(reader_->input_files()[index_]));
    }
    return *stream_;
}

FRCI::const_iterator(const FileReader* reader, int index):
    reader_(reader), index_(index)
{ }

FRCI FileReader::begin() const {
    return FRCI(this, 0);
}

FRCI FileReader::end() const {
    return FRCI(this, input_files().size());
}

bool FileReader::empty() const {
    return input_files().empty();
}

void FileReader::set_input_file(const std::string& input_file) {
    input_files_.clear();
    input_files_.push_back(input_file);
}

std::istream& FileReader::input() const {
    if (empty()) {
        throw Exception("FileReader::input() called of empty file reader");
    }
    if (!stream_) {
        stream_.reset(input_stream(input_files()[0]));
    }
    return *stream_;
}

std::istream* FileReader::input_stream(const std::string& file_name) const {
    return new std::ifstream(file_name.c_str());
}

}

