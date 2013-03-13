/*
 * bloomrepeats, Find genomic repeats, using Bloom filter based prefiltration
 * Copyright (C) 2012 Boris Nagaev
 *
 * See the LICENSE file for terms of use.
 */

#include "read_block_set.hpp"
#include "BlockSet.hpp"
#include "Block.hpp"
#include "Fragment.hpp"
#include "AlignmentRow.hpp"
#include "throw_assert.hpp"

namespace bloomrepeats {

BlockSetFastaReader::BlockSetFastaReader(BlockSet& block_set,
        std::istream& input, bool keep_alignment, RowType row_type):
    FastaReader(input), block_set_(block_set),
    keep_alignment_(keep_alignment),
    row_(0),
    row_type_(row_type)
{ }

void BlockSetFastaReader::new_sequence(const std::string& name,
                                       const std::string& description) {
    Fragment* f = block_set_.fragment_from_id(name);
    BOOST_ASSERT_MSG(f, ("No sequence or wrong position of fragment '" +
                         name + "'").c_str());
    std::string block_name = block_set_.block_from_description(description);
    BOOST_ASSERT(!block_name.empty());
    Block* block = name2block_[block_name];
    if (!block) {
        block = new Block(block_name);
        name2block_[block_name] = block;
        block_set_.insert(block);
    }
    block->insert(f);
    if (keep_alignment_) {
        row_ = AlignmentRow::new_row(row_type_);
        f->set_row(row_);
    }
}

void BlockSetFastaReader::grow_sequence(const std::string& data) {
    if (keep_alignment_) {
        BOOST_ASSERT(row_);
        row_->grow(data);
    }
}

}

