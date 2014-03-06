/*
 * bloomrepeats, Find genomic repeats, using Bloom filter based prefiltration
 * Copyright (C) 2012 Boris Nagaev
 *
 * See the LICENSE file for terms of use.
 */

#include <boost/foreach.hpp>
#include <boost/algorithm/string/split.hpp>
#include <boost/algorithm/string/classification.hpp>

#include "read_block_set.hpp"
#include "BlockSet.hpp"
#include "Block.hpp"
#include "Sequence.hpp"
#include "Fragment.hpp"
#include "AlignmentRow.hpp"
#include "complement.hpp"
#include "key_value.hpp"
#include "Exception.hpp"
#include "throw_assert.hpp"

namespace bloomrepeats {

BlockSetFastaReader::BlockSetFastaReader(BlockSet& block_set,
        std::istream& input, RowType row_type,
        SequenceType seq_type):
    FastaReader(input),
    unknown_bs_allowed_(true),
    keep_alignment_(true),
    row_type_(row_type),
    seq_type_(seq_type),
    fragment_(0),
    sequence_(0) {
    set_block_set("target", &block_set);
}

void BlockSetFastaReader::set_block_set(const std::string& name,
                                        BlockSet* block_set) {
    name2block_set_[name] = block_set;
}

void BlockSetFastaReader::new_sequence(const std::string& name,
                                       const std::string& description) {
    block_sets_.clear();
    rows_.clear();
    fragment_ = 0;
    sequence_ = 0;
    std::string block_set_name = extract_value(description, "set");
    if (block_set_name.empty()) {
        block_set_name = "target";
    }
    if (block_set_name == "all") {
        BOOST_FOREACH (const Name2BlockSet::value_type& n_bs, name2block_set_) {
            BlockSet* bs = n_bs.second;
            if (bs) {
                block_sets_.push_back(bs);
            }
        }
    } else {
        using namespace boost::algorithm;
        std::vector<std::string> names;
        split(names, block_set_name, is_any_of(","));
        BOOST_FOREACH (const std::string& block_set_name_1, names) {
            BlockSet* bs = name2block_set_[block_set_name_1];
            if (bs) {
                block_sets_.push_back(bs);
            } else if (!unknown_bs_allowed()) {
                throw Exception("Unknown block set '" + block_set_name_1 + "'");
            }
        }
    }
    if (block_sets_.empty()) {
        return;
    }
    std::string seq_name = Fragment::seq_name_from_id(name);
    if (seq_name.empty()) {
        // must be whole sequence
        SequencePtr seq = Sequence::new_sequence(seq_type_);
        seq->set_name(name);
        seq->set_description(description);
        BOOST_FOREACH (BlockSet* bs, block_sets_) {
            bs->add_sequence(seq);
        }
        sequence_ = seq.get();
        return;
    }
    SequencePtr seq = block_sets_.front()->seq_from_name(seq_name);
    if (!seq) {
        seq = Sequence::new_sequence(seq_type_);
        seq->set_name(seq_name);
        BOOST_FOREACH (BlockSet* bs, block_sets_) {
            bs->add_sequence(seq);
        }
        seqs_from_frags_.insert(seq.get());
    }
    std::string block_name = BlockSet::block_from_description(description);
    BOOST_ASSERT_MSG(!block_name.empty(), seq_name.c_str());
    keep_alignment_ = ((" " + description + " ").find("norow") ==
                       std::string::npos);
    Fragment* fragment;
    BOOST_FOREACH (BlockSet* bs, block_sets_) {
        Fragment* f = block_sets_.front()->fragment_from_id(name);
        BOOST_ASSERT_MSG(f, ("No sequence or wrong position of fragment '" +
                             name + "'").c_str());
        Block* block = name2block_[bs][block_name];
        if (!block) {
            block = new Block(block_name);
            name2block_[bs][block_name] = block;
            bs->insert(block);
        }
        block->insert(f);
        if (keep_alignment_) {
            AlignmentRow* row = AlignmentRow::new_row(row_type_);
            f->set_row(row);
            rows_.push_back(row);
        }
        fragment = f; // one of them, no matter which
    }
    if (seqs_from_frags_.find(fragment->seq()) != seqs_from_frags_.end()) {
        used_np_ = 0;
        fragment_ = fragment;
        sequence_ = fragment->seq();
    }
}

void BlockSetFastaReader::grow_sequence(const std::string& data) {
    if (block_sets_.empty()) {
        return;
    }
    if (sequence_) {
        std::string data_copy(data);
        Sequence::to_atgcn(data_copy);
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
    if (keep_alignment_) {
        BOOST_FOREACH (AlignmentRow* row, rows_) {
            row->grow(data);
        }
    }
}

}

