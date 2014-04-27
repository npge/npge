/*
 * bloomrepeats, Find genomic repeats, using Bloom filter based prefiltration
 * Copyright (C) 2014 Boris Nagaev
 *
 * See the LICENSE file for terms of use.
 */

#include <vector>
#include <map>
#include <set>
#include <algorithm>
#include <boost/foreach.hpp>

#include "SimilarAligner.hpp"
#include "throw_assert.hpp"
#include "global.hpp"

namespace bloomrepeats {

typedef std::vector<int> Ints;

struct Alignment {
    const Strings& seqs;
    Strings aligned;
    Ints pos;
    const int size;
    int mismatch_check;
    int gap_check;
    int aligned_check;

    Alignment(const Strings& s):
        seqs(s), size(s.size()) {
        aligned.resize(size);
        pos.resize(size);
    }
};

static bool equal_length(const Alignment& aln) {
    int length = aln.aligned.front().length();
    for (int i = 1; i < aln.size; i++) {
        if (aln.aligned[i].length() != length) {
            return false;
        }
    }
    return true;
}

static bool is_stop(const Alignment& aln, int shift = 0) {
    for (int i = 0; i < aln.size; i++) {
        if (aln.pos[i] + shift >= aln.seqs[i].size()) {
            return true;
        }
    }
    return false;
}

static void append_cols(Alignment& aln, int cols = 1) {
    for (int i = 0; i < aln.size; i++) {
        int& p = aln.pos[i];
        const std::string& seq = aln.seqs[i];
        std::string& a = aln.aligned[i];
        for (int j = 0; j < cols; j++) {
            a += seq[p];
            p += 1;
        }
    }
}

static void append_gaps(Alignment& aln) {
    int max_l = 0;
    for (int i = 0; i < aln.size; i++) {
        max_l = std::max(max_l, int(aln.aligned[i].length()));
    }
    for (int i = 0; i < aln.size; i++) {
        aln.aligned[i].resize(max_l, '-');
    }
    ASSERT_TRUE(equal_length(aln));
}

static void append_all(Alignment& aln) {
    for (int i = 0; i < aln.size; i++) {
        int& p = aln.pos[i];
        std::string tail = aln.seqs[i].substr(p);
        aln.aligned[i] += tail;
        p += tail.size();
    }
    append_gaps(aln);
}

static bool is_equal(const Ints& pos, const Alignment& aln,
                     int shift = 0, int cols = 1) {
    for (int j = 0; j < cols; j += 1) {
        int p = pos.front() + shift;
        char c = aln.seqs.front()[p];
        for (int i = 1; i < aln.size; i++) {
            const std::string& seq = aln.seqs[i];
            int p = pos[i] + shift;
            if (seq[p] != c) {
                return false;
            }
        }
    }
    return true;
}

static bool is_equal(const Alignment& aln,
                     int shift = 0, int cols = 1) {
    return is_equal(aln.pos, aln, shift, cols);
}

static bool is_mismatch(const Alignment& aln) {
    return !is_stop(aln, aln.mismatch_check) &&
           is_equal(aln, 1, aln.mismatch_check);
}

static bool try_mismatch(Alignment& aln) {
    if (is_mismatch(aln)) {
        append_cols(aln, aln.mismatch_check + 1);
        return true;
    } else {
        return false;
    }
}

static bool append_chars(Alignment& aln, int i, int cols = 1) {
    int& p = aln.pos[i];
    for (int j = 0; j < cols; j++) {
        aln.aligned[i] += aln.seqs[i][p];
        p += 1;
    }
}

static bool try_gap_char(Alignment& aln, char c) {
    Ints equal_pos((aln.size));
    for (int i = 0; i < aln.size; i++) {
        int p = aln.pos[i];
        bool match_this = (aln.seqs[i][p] == c);
        bool match_next = (aln.seqs[i][p + 1] == c);
        if (match_this == match_next) {
            // true, true or false,false
            return false;
        } else if (match_this) {
            equal_pos[i] = p;
        } else if (match_next) {
            equal_pos[i] = p + 1;
        }
    }
    if (is_equal(equal_pos, aln, /* shift */ 0, aln.gap_check)) {
        for (int i = 0; i < aln.size; i++) {
            if (equal_pos[i] == aln.pos[i] + 1) {
                append_chars(aln, i);
            }
        }
        append_gaps(aln);
        append_cols(aln, aln.gap_check);
        return true;
    } else {
        return false;
    }
}

static bool try_gap(Alignment& aln) {
    if (is_stop(aln, aln.gap_check)) {
        return false;
    }
    std::set<char> chars;
    for (int i = 0; i < aln.size; i++) {
        int p = aln.pos[i];
        chars.insert(aln.seqs[i][p]);
    }
    BOOST_FOREACH (char c, chars) {
        if (try_gap_char(aln, c)) {
            return true;
        } else {
            return false;
        }
    }
}

static int min_tail(const Alignment& aln) {
    int mt = aln.seqs.front().size() - aln.pos.front();
    for (int i = 1; i < aln.size; i++) {
        int t = aln.seqs[i].size() - aln.pos[i];
        mt = std::min(mt, t);
    }
    return mt;
}

typedef std::map<int, int> Seq2Pos;
typedef std::map<std::string, Seq2Pos> Found;

static std::string find_best_word(const Alignment& aln,
                                  Found& ff, int shift) {
    typedef std::set<std::string> StringSet;
    StringSet words;
    std::string best_word;
    for (int i = 0; i < aln.size; i++) {
        int p = aln.pos[i] + shift;
        std::string word = aln.seqs[i].substr(p, aln.aligned_check);
        words.insert(word);
        Seq2Pos& s2p = ff[word];
        if (s2p.find(i) == s2p.end()) {
            s2p[i] = shift;
        }
        if (s2p.size() == aln.size) {
            best_word = word;
        }
    }
    if (words.size() == 1) {
        // same word with shift
        best_word = *words.begin();
        for (int i = 0; i < aln.size; i++) {
            ff[best_word][i] = shift;
        }
    }
    return best_word;
}

static void append_aligned(Alignment& aln, const Seq2Pos& s2p) {
    for (int i = 0; i < aln.size; i++) {
        append_chars(aln, i, s2p.find(i)->second);
    }
    append_gaps(aln);
    append_cols(aln, aln.aligned_check);
}

static bool try_aligned(Alignment& aln) {
    Found ff;
    int max_shift = min_tail(aln) - aln.aligned_check;
    for (int shift = 0; shift < max_shift; shift++) {
        std::string best_word = find_best_word(aln, ff, shift);
        if (!best_word.empty()) {
            append_aligned(aln, ff[best_word]);
            return true;
        }
    }
    return false;
}

static bool pos_less(const Ints& a, const Ints& b) {
    int size = a.size();
    ASSERT_EQ(b.size(), size);
    for (int i = 0; i < size; i++) {
        if (a[i] >= b[i]) {
            return false;
        }
    }
    return true;
}

// try aligned end
static void append_end(Alignment& aln) {
    Ints end_pos((aln.size));
    for (int i = 0; i < aln.size; i++) {
        end_pos[i] = aln.seqs[i].size() - 1;
    }
    while (pos_less(aln.pos, end_pos) && is_equal(end_pos, aln)) {
        for (int i = 0; i < aln.size; i++) {
            end_pos[i] -= 1;
        }
    }
    for (int i = 0; i < aln.size; i++) {
        int cols = end_pos[i] - aln.pos[i];
        append_chars(aln, i, cols);
    }
    append_gaps(aln);
    append_all(aln);
}

static void process_cols(Alignment& aln) {
    for (int i = 0; i < aln.size; i++) {
        if (aln.seqs[i].empty()) {
            append_all(aln);
            return;
        }
    }
    while (true) {
        if (is_stop(aln)) {
            append_all(aln);
            return;
        } else if (is_equal(aln)) {
            append_cols(aln);
        } else if (try_mismatch(aln)) {
            // ok
        } else if (try_gap(aln)) {
            // ok
        } else if (try_aligned(aln)) {
            // ok
        } else {
            append_end(aln);
            return;
        }
        ASSERT_TRUE(equal_length(aln));
    }
}

void SimilarAligner::similar_aligner(Strings& seqs,
                                     int mismatch_check,
                                     int gap_check,
                                     int aligned_check) {
    if (seqs.empty()) {
        return;
    }
    Alignment aln((seqs));
    aln.mismatch_check = mismatch_check;
    aln.gap_check = gap_check;
    aln.aligned_check = aligned_check;
    process_cols(aln);
    for (int i = 0; i < aln.size; i++) {
        ASSERT_GTE(aln.aligned[i].size(), aln.seqs[i].size());
    }
    ASSERT_TRUE(equal_length(aln));
    seqs.swap(aln.aligned);
}

SimilarAligner::SimilarAligner() {
    add_opt("mismatch-check",
            "Min number of equal columns after single mismatch",
            MISMATCH_CHECK);
    add_opt("gap-check",
            "Min number of equal columns after single gap",
            GAP_CHECK);
    add_opt("aligned-check", "Min equal aligned part",
            ALIGNED_CHECK);
}

void SimilarAligner::align_seqs_impl(Strings& seqs) const {
    int mismatch_check = opt_value("mismatch-check").as<int>();
    int gap_check = opt_value("gap-check").as<int>();
    int aligned_check = opt_value("aligned-check").as<int>();
    similar_aligner(seqs, mismatch_check, gap_check, aligned_check);
}

const char* SimilarAligner::name_impl() const {
    return "Align blocks with high similarity";
}

}

