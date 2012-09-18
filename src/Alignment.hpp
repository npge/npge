/*
 * bloomrepeats, Find genomic repeats, using Bloom filter based prefiltration
 * Copyright (C) 2012 Boris Nagaev
 *
 * See the LICENSE file for terms of use.
 */

#ifndef BR_ALIGNMENT_HPP
#define BR_ALIGNMENT_HPP

#include <map>
#include <string>

#include "global.hpp"

namespace bloomrepeats {

class Alignment {
public:
    Alignment();

    int add_fragment(FragmentPtr fragment, const std::string& alignment_string);

    void remove_fragment(int index);

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
    typedef std::map<int, int> Pos2Pos;
    struct BiMap {
        Pos2Pos fragment_to_alignment;
        Pos2Pos alignment_to_fragment;
    };
    typedef std::map<int, BiMap> Maps;
    typedef std::map<FragmentPtr, int> Fragment2Index;
    typedef std::map<int, FragmentPtr> Index2Fragment;

    Maps data_;
    Fragment2Index fragment_to_index_;
    Index2Fragment index_to_fragment_;
    int length_;
};

}

#endif

