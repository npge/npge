/*
 * bloomrepeats, Find genomic repeats, using Bloom filter based prefiltration
 * Copyright (C) 2012 Boris Nagaev
 *
 * See the LICENSE file for terms of use.
 */

#include "read_block_set.hpp"
#include "BlockSet.hpp"
#include "Block.hpp"
#include "Sequence.hpp"
#include "Fragment.hpp"
#include "AlignmentRow.hpp"
#include "complement.hpp"
#include "throw_assert.hpp"

namespace bloomrepeats {

BlockSetFastaReader::BlockSetFastaReader(BlockSet& block_set,
        std::istream& input, bool keep_alignment, RowType row_type,
        SequenceType seq_type):
    FastaReader(input), block_set_(block_set),
    keep_alignment_(keep_alignment),
    row_type_(row_type),
    seq_type_(seq_type),
    row_(0),
    fragment_(0),
    sequence_(0)
{ }

void BlockSetFastaReader::new_sequence(const std::string& name,
                                       const std::string& description) {
    fragment_ = 0;
    std::string seq_name = Fragment::seq_name_from_id(name);
    if (seq_name.empty()) {
        // must be whole sequence
        SequencePtr seq = Sequence::new_sequence(seq_type_);
        seq->set_name(name);
        seq->set_description(description);
        block_set_.add_sequence(seq);
        sequence_ = seq.get();
        row_ = 0;
        return;
    }
    SequencePtr seq = block_set_.seq_from_name(seq_name);
    if (!seq) {
        seq = Sequence::new_sequence(seq_type_);
        seq->set_name(seq_name);
        block_set_.add_sequence(seq);
        seqs_from_frags_.insert(seq.get());
    }
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
    if (seqs_from_frags_.find(f->seq()) != seqs_from_frags_.end()) {
        used_np_ = 0;
        fragment_ = f;
        sequence_ = f->seq();
    }
}

void BlockSetFastaReader::grow_sequence(const std::string& data) {
    if (sequence_) {
        std::string data_copy(data);
        Sequence::to_atgc(data_copy);
        size_t min_pos;
        if (fragment_) {
            if (fragment_->ori() == 1) {
                min_pos = fragment_->min_pos() + used_np_;
            } else {
                BOOST_ASSERT(fragment_->ori() == -1);
                size_t shift = data_copy.length() + used_np_ - 1;
                min_pos = fragment_->max_pos() - shift;
                complement(data_copy);
            }
        } else {
            // must be whole sequence
            min_pos = sequence_->size();
        }
        sequence_->map_from_string(data_copy, min_pos);
        used_np_ += data_copy.length();
    }
    if (keep_alignment_ && row_) {
        row_->grow(data);
    }
}

}

