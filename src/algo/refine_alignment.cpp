/*
 * bloomrepeats, Find genomic repeats, using Bloom filter based prefiltration
 * Copyright (C) 2014 Boris Nagaev
 *
 * See the LICENSE file for terms of use.
 */

#include <algorithm>

#include "refine_alignment.hpp"
#include "throw_assert.hpp"

namespace bloomrepeats {

struct PosProps {
    bool gap;
    bool other;
    int matches;

    PosProps(const Strings& aligned, int i, int j, char c):
        gap(false), other(false), matches(0) {
        int size = aligned.size();
        for (int i1 = 0; i1 < size; i1++) {
            if (i1 != i) {
                char c1 = aligned[i1][j];
                if (c1 == '-') {
                    gap = true;
                } else if (c1 == c) {
                    matches += 1;
                } else {
                    other = true;
                }
            }
        }
    }
};

static bool can_move(const Strings& aligned, int i, int from, int to) {
    int size = aligned.size();
    int length = aligned.front().size();
    ASSERT_GTE(from, 0);
    ASSERT_LT(from, length);
    ASSERT_GTE(to, 0);
    ASSERT_LT(to, length);
    char c = aligned[i][from];
    ASSERT_NE(c, '-');
    ASSERT_EQ(aligned[i][to], '-');
    PosProps from_pos((aligned), i, from, c);
    PosProps to_pos((aligned), i, to, c);
    if (to_pos.matches == 0) {
        return false;
    }
    if (!from_pos.other) {
        return false;
    }
    if (to_pos.other && from_pos.matches) {
        return false;
    }
    return true;
}

static bool move_chars(Strings& aligned) {
    bool result = false;
    int size = aligned.size();
    int length = aligned.front().size();
    for (int i = 0; i < size; i++) {
        std::string& row = aligned[i];
        int before = -1;
        for (int j = 0; j < length; j++) {
            char c = row[j];
            if (c == '-' && j > 0 && before == -1) {
                before = j - 1;
            }
            if (c != '-' && before != -1) {
                int after = j;
                if (can_move(aligned, i, before, after - 1)) {
                    std::swap(row[before], row[after - 1]);
                    result = true;
                    before = -1;
                } else if (can_move(aligned, i, after, before + 1)) {
                    std::swap(row[after], row[before + 1]);
                    result = true;
                    before += 1;
                } else {
                    before = -1;
                }
            }
        }
    }
    return result;
}

static bool is_pure_gap(const Strings& aligned, int j) {
    int size = aligned.size();
    for (int i = 0; i < size; i++) {
        if (aligned[i][j] != '-') {
            return false;
        }
    }
    return true;
}

static void remove_pure_gaps(Strings& aligned) {
    int size = aligned.size();
    int length = aligned.front().size();
    Strings new_aligned((size));
    for (int j = 0; j < length; j++) {
        if (!is_pure_gap(aligned, j)) {
            for (int i = 0; i < size; i++) {
                new_aligned[i] += aligned[i][j];
            }
        }
    }
    aligned.swap(new_aligned);
}

void refine_alignment(Strings& aligned) {
    if (aligned.empty()) {
        return;
    }
    while (move_chars(aligned)) {
    }
    remove_pure_gaps(aligned);
}

}

