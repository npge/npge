/*
 * bloomrepeats, Find genomic repeats, using Bloom filter based prefiltration
 * Copyright (C) 2012 Boris Nagaev
 *
 * See the LICENSE file for terms of use.
 */

#include <cctype>
#include <algorithm>
#include <boost/assert.hpp>
#include <boost/foreach.hpp>

#include "Alignment.hpp"
#include "Fragment.hpp"
#include "BlockSet.hpp"
#include "Block.hpp"
#include "FastaReader.hpp"

namespace bloomrepeats {

Alignment::Alignment():
    length_(0), block_(0)
{ }

Alignment::~Alignment() {
    BOOST_FOREACH (Rows::value_type& index_and_row, rows_) {
        AlignmentRow* row = index_and_row.second;
        delete row;
    }
    rows_.clear();
    fragment_to_index_.clear();
}

int Alignment::add_row(Fragment* fragment,
                       const std::string& alignment_string) {
    int index = rows_.size();
    AlignmentRow* row = new AlignmentRow(fragment, alignment_string);
    rows_[index] = row;
    fragment_to_index_[fragment] = index;
    length_ = std::max(length_, row->length());
    return index;
}

int Alignment::add_fragment(Fragment* fragment) {
    int index = rows_.size();
    AlignmentRow* row = new AlignmentRow(fragment, "");
    rows_[index] = row;
    fragment_to_index_[fragment] = index;
    return index;
}

void Alignment::grow_row(int index, const std::string& alignment_string) {
    AlignmentRow* row = rows_[index];
    row->grow(alignment_string);
    length_ = std::max(length_, row->length());
}

void Alignment::remove_row(int index) {
    fragment_to_index_.erase(fragment_at(index));
    AlignmentRow* row = rows_[index];
    bool recalculate_length = row->length() == length();
    delete row;
    rows_.erase(index);
    if (!rows_.empty()) {
        int last_index = rows_.size();
        rows_[index] = rows_[last_index];
        rows_.erase(last_index);
    }
    if (recalculate_length) {
        length_ = 0;
        BOOST_FOREACH (Rows::value_type& index_and_row, rows_) {
            AlignmentRow* row = index_and_row.second;
            length_ = std::max(length_, int(row->length()));
        }
    }
}

int Alignment::index_of(Fragment* fragment) const {
    Fragment2Index::const_iterator it = fragment_to_index_.find(fragment);
    if (it == fragment_to_index_.end()) {
        return -1;
    } else {
        return it->second;
    }
}

Fragment* Alignment::fragment_at(int index) const {
    Rows::const_iterator it = rows_.find(index);
    if (it == rows_.end()) {
        return 0;
    } else {
        return it->second->fragment();
    }
}

int Alignment::map_to_alignment(int index, int fragment_pos) const {
    Rows::const_iterator it = rows_.find(index);
    if (it == rows_.end()) {
        return -1;
    } else {
        const AlignmentRow* row = it->second;
        return row->map_to_alignment(fragment_pos);
    }
}

int Alignment::map_to_fragment(int index, int align_pos) const {
    Rows::const_iterator it = rows_.find(index);
    if (it == rows_.end()) {
        return -1;
    } else {
        const AlignmentRow* row = it->second;
        return row->map_to_fragment(align_pos);
    }
}

int Alignment::nearest_in_fragment(int index, int align_pos) const {
    Rows::const_iterator it = rows_.find(index);
    if (it == rows_.end()) {
        return -1;
    } else {
        const AlignmentRow* row = it->second;
        return row->nearest_in_fragment(align_pos);
    }
}

int Alignment::size() const {
    return rows_.size();
}

class AlignmentFastaReader : public FastaReader {
public:
    AlignmentFastaReader(Alignment& alignment, std::istream& input):
        FastaReader(input), alignment_(alignment) {
        for (int index = 0; index < alignment.size(); index++) {
            Fragment* fragment = alignment.fragment_at(index);
            id2fragment_[fragment->id()] = fragment;
        }
        index_ = -1;
    }

    void new_sequence(const std::string& name, const std::string& description) {
        Fragment* fragment = id2fragment_[name];
        if (!fragment && alignment_.block_set()) {
            fragment = alignment_.block_set()->fragment_from_id(name);
            BOOST_ASSERT(fragment);
            id2fragment_[fragment->id()] = fragment;
            alignment_.add_fragment(fragment);
            if (!alignment_.block_) {
                Block* block = new Block;
                block->set_name(BlockSet::block_from_description(description));
                alignment_.block_set()->insert(block);
                alignment_.block_ = block;
            }
            alignment_.block_->insert(fragment);
        }
        BOOST_ASSERT(fragment);
        index_ = alignment_.index_of(fragment);
        BOOST_ASSERT(index_ != -1);
    }

    void grow_sequence(const std::string& data) {
        BOOST_ASSERT(index_ != -1);
        alignment_.grow_row(index_, data);
    }

private:
    Alignment& alignment_;
    std::map<std::string, Fragment*> id2fragment_;
    int index_;
};

std::istream& operator>>(std::istream& input, Alignment& alignment) {
    AlignmentFastaReader reader(alignment, input);
    reader.read_until_empty_line();
    return input;
}

std::ostream& operator<<(std::ostream& o, const Alignment& alignment) {
    for (int index = 0; index < alignment.size(); index++) {
        AlignmentRow* row = alignment.rows_.find(index)->second;
        o << *row;
    }
    return o;
}

}

