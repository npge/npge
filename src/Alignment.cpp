/*
 * bloomrepeats, Find genomic repeats, using Bloom filter based prefiltration
 * Copyright (C) 2012 Boris Nagaev
 *
 * See the LICENSE file for terms of use.
 */

#include <cctype>
#include <algorithm>
#include <boost/assert.hpp>

#include "Alignment.hpp"
#include "Fragment.hpp"

namespace bloomrepeats {

Alignment::Alignment():
    length_(0)
{ }

int Alignment::add_fragment(FragmentPtr fragment,
                            const std::string& alignment_string) {
    length_ = std::max(length_, int(alignment_string.length()));
    int index = data_.size();
    index_to_fragment_[index] = fragment;
    fragment_to_index_[fragment] = index;
    int non_gaps_count = 0;
    int fragment_pos = 0;
    int align_pos = 0;
    while (non_gaps_count < fragment->length()) {
        if (alignment_string[align_pos] == '-') { // FIXME other gap symbols
            align_pos += 1;
        } else {
            BOOST_ASSERT(tolower(fragment->raw_at(fragment_pos)) ==
                         tolower(alignment_string[align_pos]));
            data_[index].fragment_to_alignment[fragment_pos] = align_pos;
            data_[index].alignment_to_fragment[align_pos] = fragment_pos;
        }
    }
    return index;
}

void Alignment::remove_fragment(int index) {
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
    Index2Fragment::const_iterator it = index_to_fragment_.find(index);
    if (it == index_to_fragment_.end()) {
        return 0;
    } else {
        return it->second;
    }
}

int Alignment::map_to_alignment(int index, int fragment_pos) const {
    Maps::const_iterator it = data_.find(index);
    if (it == data_.end()) {
        return -1;
    } else {
        const Pos2Pos& fragment_to_alignment = it->second.fragment_to_alignment;
        Pos2Pos::const_iterator it2 = fragment_to_alignment.find(fragment_pos);
        if (it2 == fragment_to_alignment.end()) {
            return -1;
        } else {
            return it2->second;
        }
    }
}

int Alignment::map_to_fragment(int index, int align_pos) const {
    Maps::const_iterator it = data_.find(index);
    if (it == data_.end()) {
        return -1;
    } else {
        const Pos2Pos& alignment_to_fragment = it->second.alignment_to_fragment;
        Pos2Pos::const_iterator it2 = alignment_to_fragment.find(align_pos);
        if (it2 == alignment_to_fragment.end()) {
            return -1;
        } else {
            return it2->second;
        }
    }
}

int Alignment::nearest_in_fragment(int index, int align_pos) const {
    // FIXME do smth with this
    for (int distance = 0; distance < length(); distance++) {
        for (int ori = -1; ori <= 1; ori += 2) {
            int new_align_pos = align_pos + ori * distance;
            if (map_to_fragment(index, new_align_pos) != -1) {
                return map_to_fragment(index, new_align_pos);
            }
        }
    }
    return -1;
}

int Alignment::size() const {
    return data_.size();
}

}

