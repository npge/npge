/*
 * NPG-explorer, Nucleotide PanGenome explorer
 * Copyright (C) 2012-2015 Boris Nagaev
 *
 * See the LICENSE file for terms of use.
 */

#include <vector>

#include "FragmentDistance.hpp"
#include "Block.hpp"
#include "Fragment.hpp"
#include "AlignmentRow.hpp"
#include "Exception.hpp"

namespace npge {

static int substitution(char a, char b) {
    // TODO move to util/, see PairAligner
    return (a == b && a != 'N') ? 0 : 1;
}

static void equal_at(const Fragment* a, const Fragment* b,
                     pos_t pos, bool& eq, bool& gap) {
    char ac = a->alignment_at(pos);
    char bc = b->alignment_at(pos);
    eq = (ac == bc);
    gap = (ac == '\0');
}

FragmentDistance::Distance FragmentDistance::fragment_distance(
    const Fragment* a, const Fragment* b) const {
    TimeIncrementer ti(this);
    AlignmentRow* ar = a->row();
    AlignmentRow* br = b->row();
    if (!ar || !br) {
        throw Exception("Fragment without alignment");
    }
    if (ar->length() != br->length()) {
        throw Exception("Alignment rows of different lengths");
    }
    int length = ar->length();
    Distance result;
    result.total = length;
    result.penalty = 0;
    if (length < 3) {
        return result;
    }
    // cycle buffer
    bool eq[3];
    bool gap[3];
    for (int i = 0; i < 2; i++) {
        equal_at(a, b, i, eq[i], gap[i]);
    }
    for (int i = 2; i < length; i++) {
        // i is index of third letter (added)
        int prev = (i - 2) % 3;
        int curr = (i - 1) % 3;
        int next = i % 3;
        equal_at(a, b, i, eq[next], gap[next]);
        bool prev_good = eq[prev] && !gap[prev];
        bool curr_good = !eq[curr];
        bool next_good = eq[next] && !gap[next];
        if (prev_good && curr_good && next_good) {
            result.penalty += 1;
        }
    }
    return result;
}

double FragmentDistance::Distance::ratio() const {
    return double(penalty) / double(total);
}

FragmentDistance::FragmentDistance() {
    set_opt_prefix("distance-");
    declare_bs("target", "Target blockset");
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

const char* FragmentDistance::name_impl() const {
    return "Calculate distance between fragments in block";
}

}

