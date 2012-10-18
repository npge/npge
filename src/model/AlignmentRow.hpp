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

// TODO memory-friendly implementation

class AlignmentRow {
public:
    AlignmentRow(Fragment* fragment);

    virtual void grow(const std::string& alignment_string);

    virtual void bind(int fragment_pos, int align_pos) = 0;

    virtual int map_to_alignment(int fragment_pos) const = 0;

    virtual int map_to_fragment(int align_pos) const = 0;

    int length() const {
        return length_;
    }

    Fragment* fragment() const {
        return fragment_;
    }

    virtual int nearest_in_fragment(int align_pos) const;

    virtual void print_alignment_string(std::ostream& o) const;

protected:
    void set_length(int length) {
        length_ = length;
    }

private:
    int length_;
    Fragment* fragment_;
};

class MapAlignmentRow : public AlignmentRow {
public:
    MapAlignmentRow(Fragment* fragment, const std::string& alignment_string);

    void bind(int fragment_pos, int align_pos);

    int map_to_alignment(int fragment_pos) const;

    int map_to_fragment(int align_pos) const;

private:
    typedef std::map<int, int> Pos2Pos;

    Pos2Pos fragment_to_alignment_;
    Pos2Pos alignment_to_fragment_;
};

/** Streaming operator */
std::ostream& operator<<(std::ostream& o, const AlignmentRow& row);

}

#endif

