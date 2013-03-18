/*
 * bloomrepeats, Find genomic repeats, using Bloom filter based prefiltration
 * Copyright (C) 2012 Boris Nagaev
 *
 * See the LICENSE file for terms of use.
 */

#include <boost/foreach.hpp>

#include "block_stat.hpp"
#include "Block.hpp"
#include "Fragment.hpp"
#include "boundaries.hpp"

namespace bloomrepeats {

AlignmentStat::AlignmentStat():
    ident_nogap(0),
    ident_gap(0),
    noident_nogap(0),
    noident_gap(0),
    pure_gap(0),
    total(0),
    spreading(0)
{ }

// TODO rename Boundaries to smth
typedef Boundaries Integers;

void make_stat(AlignmentStat& stat, const Block* block) {
    stat.total = block->alignment_length();
    for (size_t pos = 0; pos < stat.total; pos++) {
        char seen_letter = 0;
        bool ident = true;
        bool gap = false;
        BOOST_FOREACH (Fragment* f, *block) {
            char c = f->alignment_at(pos);
            if (c == 0) {
                gap = true;
            } else if (seen_letter == 0) {
                seen_letter = c;
            } else if (c != seen_letter) {
                ident = false;
            }
        }
        if (seen_letter) {
            if (ident && !gap) {
                stat.ident_nogap += 1;
            } else if (ident && gap) {
                stat.ident_gap += 1;
            } else if (!ident && !gap) {
                stat.noident_nogap += 1;
            } else if (!ident && gap) {
                stat.noident_gap += 1;
            }
        } else {
            stat.pure_gap += 1;
        }
    }
    Integers lengths;
    BOOST_FOREACH (Fragment* f, *block) {
        lengths.push_back(f->length());
    }
    if (!lengths.empty()) {
        int max_length = *std::max_element(lengths.begin(), lengths.end());
        int min_length = *std::min_element(lengths.begin(), lengths.end());
        int avg_length = avg_element(lengths);
        if (avg_length == 0) {
            stat.spreading = 0;
        } else {
            stat.spreading = float(max_length - min_length) / avg_length;
        }
    }
}

float block_identity(const AlignmentStat& stat, bool allow_gaps) {
    int accepted = stat.ident_nogap;
    int total = stat.ident_nogap + stat.noident_nogap;
    if (allow_gaps) {
        accepted += stat.ident_gap;
        total += stat.stat.ident_gap + stat.noident_gap;
    }
    return total ? float(accepted) / float(total) : 0;
}

}

