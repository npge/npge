/*
 * bloomrepeats, Find genomic repeats, using Bloom filter based prefiltration
 * Copyright (C) 2012 Boris Nagaev
 *
 * See the LICENSE file for terms of use.
 */

#include <cstdlib>
#include <boost/format.hpp>
#include <boost/foreach.hpp>

#include "ExternalAligner.hpp"
#include "AlignmentRow.hpp"
#include "FastaReader.hpp"
#include "Block.hpp"
#include "Fragment.hpp"
#include "temp_file.hpp"
#include "name_to_stream.hpp"
#include "throw_assert.hpp"

namespace bloomrepeats {

class AlignmentReader : public FastaReader {
public:
    AlignmentReader(Block& block, std::istream& input, RowType type):
        FastaReader(input),
        block_(block),
        row_(0),
        row_type_(type)
    { }

    void new_sequence(const std::string& name, const std::string& description) {
        row_ = 0;
        BOOST_FOREACH (Fragment* f, block_) {
            if (f->id() == name) {
                row_ = AlignmentRow::new_row(row_type_);
                f->set_row(row_);
            }
        }
        BOOST_ASSERT(row_);
    }

    void grow_sequence(const std::string& data) {
        BOOST_ASSERT(row_);
        row_->grow(data);
    }

private:
    Block& block_;
    AlignmentRow* row_;
    RowType row_type_;
};

ExternalAligner::ExternalAligner(const std::string& cmd):
    cmd_(cmd)
{ }

void ExternalAligner::align_block(Block* block) const {
    if (block->size() <= 1) {
        BOOST_FOREACH (Fragment* f, *block) {
            AlignmentRow* row = AlignmentRow::new_row(COMPACT_ROW);
            int length = f->length();
            row->set_length(length);
            for (int i = 0; i < length; i++) {
                row->bind(i, i);
            }
            f->set_row(row);
        }
        return;
    }
    std::string input = temp_file();
    BOOST_ASSERT(!input.empty());
    std::string output = temp_file();
    BOOST_ASSERT(!output.empty());
    {
        boost::shared_ptr<std::ostream> unaligned = name_to_ostream(input);
        *unaligned << *block;
    }
    std::string cmd_string = str(boost::format(cmd()) % input % output);
    system(cmd_string.c_str());
    {
        boost::shared_ptr<std::istream> aligned = name_to_istream(output);
        // FIXME row type
        AlignmentReader reader(*block, *aligned, COMPACT_ROW);
        reader.read_all_sequences();
    }
    remove_file(input);
    remove_file(output);
}

void ExternalAligner::add_options_impl(po::options_description& desc) const {
    add_unique_options(desc)
    ("aligner-cmd", po::value<std::string>()->default_value(cmd()),
     "Template of command for external aligner")
   ;
}

void ExternalAligner::apply_options_impl(const po::variables_map& vm) {
    if (vm.count("aligner-cmd")) {
        set_cmd(vm["aligner-cmd"].as<std::string>());
    }
}

bool ExternalAligner::apply_to_block_impl(Block* block) const {
    align_block(block);
    return true;
}

const char* ExternalAligner::name_impl() const {
    return "External aligner";
}

}

