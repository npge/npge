/*
 * NPG-explorer, Nucleotide PanGenome explorer
 * Copyright (C) 2012-2015 Boris Nagaev
 *
 * See the LICENSE file for terms of use.
 */

#include <algorithm>

#include "refine_alignment.hpp"
#include "throw_assert.hpp"

namespace npge {

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
    const std::string& row = aligned[i];
    char c = row[from];
    char to_c = row[to];
    if ((c == '-') == (to_c == '-')) {
        return false;
    }
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

static bool try_move(Strings& aligned, int i, int from, int to) {
    bool result = can_move(aligned, i, from, to);
    if (result) {
        std::string& row = aligned[i];
        std::swap(row[from], row[to]);
    }
    return result;
}

// TODO DRY
static bool is_equal(const Strings& aligned, int j) {
    int size = aligned.size();
    char c = aligned.front()[j];
    for (int i = 1; i < size; i++) {
        if (aligned[i][j] != c) {
            return false;
        }
    }
    return true;
}

static bool check_movable(Strings& aligned,
                          int i, int first, int last) {
    int l = aligned.front().size();
    const std::string& row = aligned[i];
    ASSERT_EQ(row[first], row[last]);
    ASSERT_TRUE(first == 0 || row[first - 1] != row[first]);
    ASSERT_TRUE(last == l - 1 || row[last + 1] != row[last]);
    if (row[first] == '-') {
        if (first > 0 && try_move(aligned, i, first - 1, last)) {
            // aaa----bbbbb
            // aaaB----bbbb
            return true;
        }
        if (last < l - 1 && try_move(aligned, i, last + 1, first)) {
            // aaa----Abbbb
            // aaaa----bbbb
            return true;
        }
    } else {
        if (last < l - 1 && try_move(aligned, i, first, last + 1)) {
            // aaaBBBB-
            // aaacbbbb
            return true;
        }
        if (first > 0 && try_move(aligned, i, last, first - 1)) {
            // aaa-BBBB
            // aaabbbbc
            return true;
        }
        for (int j = first + 1; j <= last - 1; j++) {
            if (!is_equal(aligned, j)) {
                if (last < l - 1 &&
                        try_move(aligned, i, j, last + 1)) {
                    return true;
                }
                if (first > 0 &&
                        try_move(aligned, i, j, first - 1)) {
                    return true;
                }
            }
        }
    }
    return false;
}

static bool move_chars(Strings& aligned) {
    bool result = false;
    int size = aligned.size();
    int length = aligned.front().size();
    for (int i = 0; i < size; i++) {
        std::string& row = aligned[i];
        char repeated = row[0];
        int first = 0, last = 0;
        for (int j = 1; j < length; j++) {
            char c = row[j];
            if (c == repeated) {
                last = j;
            } else {
                result |= check_movable(aligned, i, first, last);
                repeated = row[j];
                first = j;
                last = j;
                while (first > 0 && row[first - 1] == repeated) {
                    first -= 1;
                }
            }
        }
        result |= check_movable(aligned, i, first, last);
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
        remove_pure_gaps(aligned);
    }
    remove_pure_gaps(aligned);
}

}

