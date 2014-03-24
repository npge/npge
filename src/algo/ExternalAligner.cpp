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
    if (!alignment_needed(block)) {
        return;
    }
    std::string input = tmp_file();
    BOOST_ASSERT(!input.empty());
    std::string output = tmp_file();
    BOOST_ASSERT(!output.empty());
    std::vector<Fragment*> fragments((block->begin()), block->end());
    {
        boost::shared_ptr<std::ostream> file = name_to_ostream(input);
        BOOST_FOREACH (Fragment* f, fragments) {
            *file << *f;
        }
    }
    align_file(input, output);
    std::vector<std::string> rows;
    read_alignment(rows, output);
    BOOST_ASSERT(rows.size() == fragments.size());
    for (int i = 0; i < fragments.size(); i++) {
        new CompactAlignmentRow(rows[i], fragments[i]);
    }
    remove_file(input);
    remove_file(output);
}

bool ExternalAligner::alignment_needed(Block* block) const {
    if (block->size() == 0) {
        return false;
    } else if (block->size() == 1) {
        Fragment* f = block->front();
        if (f->row() && f->row()->length() == f->length()) {
            return false;
        }
        AlignmentRow* row = AlignmentRow::new_row(COMPACT_ROW);
        int length = f->length();
        row->set_length(length);
        for (int i = 0; i < length; i++) {
            row->bind(i, i);
        }
        f->set_row(row);
        return false;
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
            return false;
        }
    }
    return true;
}

void ExternalAligner::align_file(const std::string& input,
                                 const std::string& output) const {
    std::string cmd = opt_value("aligner-cmd").as<std::string>();
    std::string cmd_string = str(boost::format(cmd) % input % output);
    int r = system(cmd_string.c_str());
    if (r) {
        throw Exception("external aligner failed with code " +
                        TO_S(r));
    }
}

class AlignmentReader : public FastaReader {
public:
    AlignmentReader(std::vector<std::string>& rows,
                    std::istream& input):
        rows_(rows),
        FastaReader(input)
    { }

    void new_sequence(const std::string& name,
                      const std::string& description) {
        rows_.push_back("");
    }

    void grow_sequence(const std::string& data) {
        BOOST_ASSERT(!rows_.empty());
        rows_.back() += data;
    }

    std::vector<std::string>& rows_;
};

void ExternalAligner::read_alignment(std::vector<std::string>& rows,
                                     const std::string& file) const {
    boost::shared_ptr<std::istream> aligned = name_to_istream(file);
    AlignmentReader reader(rows, *aligned);
    reader.read_all_sequences();
}

void ExternalAligner::process_block_impl(Block* block, ThreadData*) const {
    align_block(block);
}

const char* ExternalAligner::name_impl() const {
    return "External aligner";
}

}

