/*
 * bloomrepeats, Find genomic repeats, using Bloom filter based prefiltration
 * Copyright (C) 2012 Boris Nagaev
 *
 * See the LICENSE file for terms of use.
 */

#ifndef BR_ALIGNMENT_ROW_HPP_
#define BR_ALIGNMENT_ROW_HPP_

#include <map>
#include <string>

#include "global.hpp"

namespace bloomrepeats {

class AlignmentRow {
public:
    AlignmentRow(Fragment* fragment, const std::string& alignment_string);

    void grow(const std::string& alignment_string);

    int map_to_alignment(int fragment_pos) const;

    int map_to_fragment(int align_pos) const;

    int nearest_in_fragment(int align_pos) const;

    // TODO
    //void bind(int index, int fragment_pos, int align_pos);
    // read from and write to stream

    int length() const {
        return length_;
    }

    Fragment* fragment() const {
        return fragment_;
    }

    void print_alignment_string(std::ostream& o) const;

private:
    // TODO memory-friendly implementation
    typedef std::map<int, int> Pos2Pos;

    Pos2Pos fragment_to_alignment;
    Pos2Pos alignment_to_fragment;
    int length_;
    Fragment* fragment_;
};

/** Streaming operator */
std::ostream& operator<<(std::ostream& o, const AlignmentRow& row);

}

#endif

