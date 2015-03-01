/*
 * NPG-explorer, Nucleotide PanGenome explorer
 * Copyright (C) 2012-2015 Boris Nagaev
 *
 * See the LICENSE file for terms of use.
 */

#include "FileWriter.hpp"
#include "Processor.hpp"
#include "name_to_stream.hpp"
#include "throw_assert.hpp"

namespace npge {

FileWriter::FileWriter(Processor* processor,
                       const std::string& opt,
                       const std::string& descr,
                       bool required,
                       const std::string& dflt):
    processor_(processor), opt_(opt) {
    processor_->add_opt(opt_, descr, dflt, required);
    processor_->add_opt("remove-after", "remove file " + opt_,
                        false);
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
    reset();
}

void FileWriter::set_rand_name(bool remove_prev) {
    ASSERT_TRUE(processor_);
    set_output_file(processor_->tmp_file(), remove_prev);
}

void FileWriter::reset() {
    output_.reset();
}

void FileWriter::set_remove_after(bool value) {
    processor_->set_opt_value("remove-after", value);
}

bool FileWriter::get_remove_after() const {
    return processor_->opt_value("remove-after").as<bool>();
}

std::ostream& FileWriter::output() const {
    if (!output_ || (output_file()[0] != ':' &&
                     !file_exists(output_file()))) {
        output_ = name_to_ostream(output_file());
    }
    return *output_;
}

}

