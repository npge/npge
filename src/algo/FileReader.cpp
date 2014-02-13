/*
 * bloomrepeats, Find genomic repeats, using Bloom filter based prefiltration
 * Copyright (C) 2012 Boris Nagaev
 *
 * See the LICENSE file for terms of use.
 */

#include "FileReader.hpp"
#include "Processor.hpp"
#include "Exception.hpp"
#include "name_to_stream.hpp"
#include "boost-assert.hpp"

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
        stream_ = name_to_istream(reader_->input_files()[index_]);
    }
    return *stream_;
}

FRCI::const_iterator(const FileReader* reader, int index):
    reader_(reader), index_(index)
{ }

FileReader::FileReader(Processor* processor, const std::string& opt,
        const std::string& descr):
    processor_(processor), opt_(opt) {
    processor_->add_opt(opt_, descr, Files(), true);
}

FRCI FileReader::begin() const {
    return FRCI(this, 0);
}

FRCI FileReader::end() const {
    return FRCI(this, input_files().size());
}

bool FileReader::empty() const {
    return input_files().empty();
}

FileReader::Files FileReader::input_files() const {
    BOOST_ASSERT(processor_);
    return processor_->opt_value(opt_).as<Files>();
}

void FileReader::set_input_files(const FileReader::Files& input_files) {
    BOOST_ASSERT(processor_);
    processor_->set_opt_value(opt_, input_files);
}

void FileReader::set_input_file(const std::string& input_file) {
    Files files;
    files.push_back(input_file);
    set_input_files(files);
}

std::istream& FileReader::input() const {
    if (empty()) {
        throw Exception("FileReader::input() called of empty file reader");
    }
    if (!stream_) {
        stream_ = name_to_istream(input_files()[0]);
    }
    return *stream_;
}

}

