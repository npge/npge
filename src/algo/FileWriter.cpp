/*
 * bloomrepeats, Find genomic repeats, using Bloom filter based prefiltration
 * Copyright (C) 2012 Boris Nagaev
 *
 * See the LICENSE file for terms of use.
 */

#include <boost/filesystem.hpp>

#include "FileWriter.hpp"
#include "Processor.hpp"
#include "name_to_stream.hpp"
#include "temp_file.hpp"

namespace bloomrepeats {

FileWriter::FileWriter(Processor* processor, const std::string& opt,
        const std::string& descr, bool required):
    processor_(processor), opt_(opt) {
    processor_->add_opt(opt_, descr, std::string(), required);
    processor_->add_opt("remove-after", "remove file " + opt_, false);
}

FileWriter::~FileWriter() {
    if (get_remove_after()) {
        remove_file(output_file());
    }
}

std::string FileWriter::output_file() const {
    return processor_->opt_value(opt_).as<std::string>();
}

void FileWriter::set_output_file(const std::string& output_file,
                                 bool remove_prev) {
    if (get_remove_after() || remove_prev) {
        remove_file(this->output_file());
    }
    processor_->set_opt_value(opt_, output_file);
    output_.reset();
}

void FileWriter::set_rand_name(bool remove_prev) {
    set_output_file(temp_file(), remove_prev);
}

void FileWriter::set_remove_after(bool value) {
    processor_->set_opt_value("remove-after", value);
}

bool FileWriter::get_remove_after() const {
    return processor_->opt_value("remove-after").as<bool>();
}

std::ostream& FileWriter::output() const {
    using namespace boost::filesystem;
    if (!output_ || (output_file()[0] != ':' && !exists(output_file()))) {
        output_ = name_to_ostream(output_file());
    }
    return *output_;
}

}

