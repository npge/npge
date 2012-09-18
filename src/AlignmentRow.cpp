/*
 * bloomrepeats, Find genomic repeats, using Bloom filter based prefiltration
 * Copyright (C) 2012 Boris Nagaev
 *
 * See the LICENSE file for terms of use.
 */

#include <cctype>
#include <boost/assert.hpp>

#include "AlignmentRow.hpp"
#include "Fragment.hpp"

namespace bloomrepeats {

AlignmentRow::AlignmentRow(FragmentPtr fragment,
                           const std::string& alignment_string):
    length_(alignment_string.length()), fragment_(fragment) {
    int fragment_pos = 0;
    for (int align_pos = 0; align_pos < alignment_string.size(); align_pos++) {
        if (isalpha(alignment_string[align_pos])) {
            BOOST_ASSERT(tolower(fragment->raw_at(fragment_pos)) ==
                         tolower(alignment_string[align_pos]));
            fragment_to_alignment[fragment_pos] = align_pos;
            alignment_to_fragment[align_pos] = fragment_pos;
            fragment_pos += 1;
        }
    }
}

int AlignmentRow::map_to_alignment(int fragment_pos) const {
    Pos2Pos::const_iterator it2 = fragment_to_alignment.find(fragment_pos);
    if (it2 == fragment_to_alignment.end()) {
        return -1;
    } else {
        return it2->second;
    }
}

int AlignmentRow::map_to_fragment(int align_pos) const {
    Pos2Pos::const_iterator it2 = alignment_to_fragment.find(align_pos);
    if (it2 == alignment_to_fragment.end()) {
        return -1;
    } else {
        return it2->second;
    }
}

int AlignmentRow::nearest_in_fragment(int align_pos) const {
    // FIXME do smth with this
    for (int distance = 0; distance < length(); distance++) {
        for (int ori = -1; ori <= 1; ori += 2) {
            int new_align_pos = align_pos + ori * distance;
            if (map_to_fragment(new_align_pos) != -1) {
                return map_to_fragment(new_align_pos);
            }
        }
    }
    return -1;
}

}

