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

namespace bloomrepeats {

Alignment::Alignment():
    length_(0)
{ }

Alignment::~Alignment() {
    BOOST_FOREACH (Rows::value_type& index_and_row, data_) {
        AlignmentRow* row = index_and_row.second;
        delete row;
    }
    data_.clear();
    fragment_to_index_.clear();
}

int Alignment::add_fragment(FragmentPtr fragment,
                            const std::string& alignment_string) {
    length_ = std::max(length_, int(alignment_string.length()));
    int index = data_.size();
    data_[index] = new AlignmentRow(fragment, alignment_string);
    fragment_to_index_[fragment] = index;
    return index;
}

void Alignment::remove_fragment(int index) {
    delete data_[index];
    data_.erase(index);
    if (!data_.empty()) {
        int last_index = data_.size();
        data_[index] = data_[last_index];
        data_.erase(last_index);
    }
    // TODO change length_ if this was the longest fragment
}

int Alignment::index_of(FragmentPtr fragment) const {
    Fragment2Index::const_iterator it = fragment_to_index_.find(fragment);
    if (it == fragment_to_index_.end()) {
        return -1;
    } else {
        return it->second;
    }
}

FragmentPtr Alignment::fragment_at(int index) const {
    Rows::const_iterator it = data_.find(index);
    if (it == data_.end()) {
        return 0;
    } else {
        return it->second->fragment();
    }
}

int Alignment::map_to_alignment(int index, int fragment_pos) const {
    Rows::const_iterator it = data_.find(index);
    if (it == data_.end()) {
        return -1;
    } else {
        const AlignmentRow* row = it->second;
        return row->map_to_alignment(fragment_pos);
    }
}

int Alignment::map_to_fragment(int index, int align_pos) const {
    Rows::const_iterator it = data_.find(index);
    if (it == data_.end()) {
        return -1;
    } else {
        const AlignmentRow* row = it->second;
        return row->map_to_fragment(align_pos);
    }
}

int Alignment::nearest_in_fragment(int index, int align_pos) const {
    Rows::const_iterator it = data_.find(index);
    if (it == data_.end()) {
        return -1;
    } else {
        const AlignmentRow* row = it->second;
        return row->nearest_in_fragment(align_pos);
    }
}

int Alignment::size() const {
    return data_.size();
}

}

