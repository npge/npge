/*
 * bloomrepeats, Find genomic repeats, using Bloom filter based prefiltration
 * Copyright (C) 2012 Boris Nagaev
 *
 * See the LICENSE file for terms of use.
 */

#ifndef BR_ALIGNMENT_HPP
#define BR_ALIGNMENT_HPP

#include <iosfwd>
#include <map>
#include <string>

#include "global.hpp"
#include "AlignmentRow.hpp"

namespace bloomrepeats {

class Alignment {
public:
    Alignment();

    ~Alignment();

    int add_row(FragmentPtr fragment, const std::string& alignment_string);

    int add_fragment(FragmentPtr fragment); // with empty body

    void grow_row(int index, const std::string& alignment_string);

    void remove_row(int index);

    int index_of(FragmentPtr fragment) const;

    FragmentPtr fragment_at(int index) const;

    int map_to_alignment(int index, int fragment_pos) const;

    int map_to_fragment(int index, int align_pos) const;

    int nearest_in_fragment(int index, int align_pos) const;

    // TODO
    //void bind(int index, int fragment_pos, int align_pos);
    // read from and write to stream

    int size() const;

    int length() const {
        return length_;
    }

private:
    // TODO memory-friendly implementation
    typedef std::map<int, AlignmentRow*> Rows;
    typedef std::map<FragmentPtr, int> Fragment2Index;

    Rows rows_;
    Fragment2Index fragment_to_index_;
    int length_;

    friend std::ostream& operator<<(std::ostream&, const Alignment&);
};

/** Streaming operator.
\note Fragment list must be pre-added using Alignment::add_fragment().
    Names of sequences must correspond Fragment::id().
*/
std::istream& operator>>(std::istream& i, Alignment& alignment);

/** Streaming operator */
std::ostream& operator<<(std::ostream& o, const Alignment& alignment);

}

#endif

