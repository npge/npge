/*
 * bloomrepeats, Find genomic repeats, using Bloom filter based prefiltration
 * Copyright (C) 2014 Boris Nagaev
 *
 * See the LICENSE file for terms of use.
 */

#include <vector>
#include "FragmentDistance.hpp"
#include "Block.hpp"
#include "Fragment.hpp"
#include "AlignmentRow.hpp"
#include "Exception.hpp"

namespace bloomrepeats {

static int substitution(char a, char b) {
    // TODO move to util/, see PairAligner
    return (a == b && a != 'N') ? 0 : 1;
}

FragmentDistance::Distance FragmentDistance::fragment_distance(
    const Fragment* a, const Fragment* b) const {
    AlignmentRow* ar = a->row();
    AlignmentRow* br = b->row();
    if (!ar || !br) {
        throw Exception("Fragment without alignment");
    }
    if (ar->length() != br->length()) {
        throw Exception("Alignment rows of different lengths");
    }
    int length = ar->length();
    int distance = 0;
    int total = 0;
    bool a_gaps = false;
    bool b_gaps = false;
    for (int i = 0; i < length; i++) {
        int a_fr_pos = ar->map_to_fragment(i);
        int b_fr_pos = br->map_to_fragment(i);
        if (a_fr_pos == -1 && b_fr_pos == -1) {
            continue;
        }
        total += 1;
        if (a_fr_pos == -1) {
            if (!a_gaps) {
                a_gaps = true;
                distance += 1;
            }
        } else {
            a_gaps = false;
        }
        if (b_fr_pos == -1) {
            if (!b_gaps) {
                b_gaps = true;
                distance += 1;
            }
        } else {
            b_gaps = false;
        }
        if (a_fr_pos != -1 && b_fr_pos != -1) {
            char a_char = a->raw_at(a_fr_pos);
            char b_char = b->raw_at(b_fr_pos);
            distance += substitution(a_char, b_char);
        }
    }
    Distance result;
    result.total = total;
    result.penalty = distance;
    return result;
}

double FragmentDistance::Distance::ratio() const {
    return double(penalty) / double(total);
}

FragmentDistance::FragmentDistance() {
    set_opt_prefix("distance-");
}

void FragmentDistance::print_block(std::ostream& o, Block* block) const {
    std::vector<const Fragment*> fragments(block->begin(), block->end());
    for (int i = 0; i < fragments.size(); i++) {
        const Fragment* f1 = fragments[i];
        for (int j = i + 1; j < fragments.size(); j++) {
            const Fragment* f2 = fragments[j];
            o << block->name() << '\t';
            o << f1->id() << '\t';
            o << f2->id() << '\t';
            o << fragment_distance(f1, f2).ratio() << '\n';
        }
    }
}

void FragmentDistance::print_header(std::ostream& o) const {
    o << "block" << '\t' << "fr.1" << '\t' << "fr.2"
      << '\t' << "distance" << '\n';
}

}

