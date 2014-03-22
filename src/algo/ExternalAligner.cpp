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
#include "name_to_stream.hpp"
#include "throw_assert.hpp"
#include "Exception.hpp"
#include "to_s.hpp"

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

ExternalAligner::ExternalAligner(const std::string& cmd) {
    add_opt("aligner-cmd",
            "Template of command for external aligner", cmd);
}

struct BlockSquareLess {
    bool operator()(Block* a, Block* b) const {
        int a_length = a->front() ? a->front()->length() : 0;
        int b_length = b->front() ? b->front()->length() : 0;
        return b_length * b->size() < a_length * a->size();
    }
};

void ExternalAligner::change_blocks_impl(std::vector<Block*>& blocks) const {
    std::sort(blocks.begin(), blocks.end(), BlockSquareLess());
}

void ExternalAligner::align_block(Block* block) const {
    if (block->size() == 0) {
        return;
    } else if (block->size() == 1) {
        Fragment* f = block->front();
        if (f->row() && f->row()->length() == f->length()) {
            return;
        }
        AlignmentRow* row = AlignmentRow::new_row(COMPACT_ROW);
        int length = f->length();
        row->set_length(length);
        for (int i = 0; i < length; i++) {
            row->bind(i, i);
        }
        f->set_row(row);
        return;
    }
    if (block->front()->row()) {
        int row_length = block->front()->row()->length();
        bool all_rows = true;
        BOOST_FOREACH (Fragment* f, *block) {
            if (!f->row() || f->row()->length() != row_length) {
                all_rows = false;
                break;
            }
        }
        if (all_rows) {
            // all fragments have rows and lengthes are equal
            return;
        }
    }
    std::string input = tmp_file();
    BOOST_ASSERT(!input.empty());
    std::string output = tmp_file();
    BOOST_ASSERT(!output.empty());
    {
        boost::shared_ptr<std::ostream> unaligned = name_to_ostream(input);
        *unaligned << *block;
    }
    std::string cmd = opt_value("aligner-cmd").as<std::string>();
    std::string cmd_string = str(boost::format(cmd) % input % output);
    int r = system(cmd_string.c_str());
    if (r) {
        throw Exception("external aligner failed with code " + TO_S(r));
    }
    {
        boost::shared_ptr<std::istream> aligned = name_to_istream(output);
        // FIXME row type
        AlignmentReader reader(*block, *aligned, COMPACT_ROW);
        reader.read_all_sequences();
    }
    remove_file(input);
    remove_file(output);
}

void ExternalAligner::process_block_impl(Block* block, ThreadData*) const {
    align_block(block);
}

const char* ExternalAligner::name_impl() const {
    return "External aligner";
}

}

