/*
 * NPG-explorer, Nucleotide PanGenome explorer
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
#include "FindLowSimilar.hpp"
#include "throw_assert.hpp"
#include "global.hpp"

namespace npge {

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

static void process_cols(Alignment& aln);
static void process_seqs(Strings& seqs,
                         int mismatch_check,
                         int gap_check,
                         int aligned_check);

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
        int p = pos.front() + shift + j;
        char c = aln.seqs.front()[p];
        for (int i = 1; i < aln.size; i++) {
            const std::string& seq = aln.seqs[i];
            int p = pos[i] + shift + j;
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
        }
    }
    return false;
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
    Strings tmp_seqs((aln.size));
    for (int i = 0; i < aln.size; i++) {
        int p = aln.pos[i];
        int length = s2p.find(i)->second;
        std::string& tmp_seq = tmp_seqs[i];
        tmp_seq = aln.seqs[i].substr(p, length);
        std::reverse(tmp_seq.begin(), tmp_seq.end());
    }
    process_seqs(tmp_seqs, aln.mismatch_check,
                 aln.gap_check, aln.aligned_check);
    for (int i = 0; i < aln.size; i++) {
        std::string& tmp_row = tmp_seqs[i];
        std::reverse(tmp_row.begin(), tmp_row.end());
        aln.aligned[i] += tmp_row;
        int& p = aln.pos[i];
        int length = s2p.find(i)->second;
        p += length;
    }
}

static bool try_aligned(Alignment& aln) {
    Found ff;
    int max_shift = min_tail(aln) - aln.aligned_check;
    for (int shift = 0; shift < max_shift; shift++) {
        std::string best_word = find_best_word(aln, ff, shift);
        if (!best_word.empty()) {
            append_aligned(aln, ff[best_word]);
            append_cols(aln, aln.aligned_check);
            return true;
        }
    }
    return false;
}

static bool pos_less(const Ints& a, const Ints& b, int shift = 0) {
    int size = a.size();
    ASSERT_EQ(b.size(), size);
    for (int i = 0; i < size; i++) {
        if (a[i] >= b[i] + shift) {
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
    while ((pos_less(aln.pos, end_pos) &&
            is_equal(end_pos, aln)) ||
            (pos_less(aln.pos, end_pos, -1) &&
             is_equal(end_pos, aln, -1))) {
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

typedef FindLowSimilar::Region Region;
typedef std::vector<Region> Regions;

static void append_region(Strings& aligned, const Alignment& aln,
                          const Region& r) {
    for (int i = 0; i < aln.size; i++) {
        aligned[i] += aln.seqs[i].substr(r.start_, r.length());
    }
}

static void filter_out_gaps(Strings& aligned) {
    int size = aligned.size();
    for (int i = 0; i < size; i++) {
        std::string& a = aligned[i];
        a.erase(std::remove(a.begin(), a.end(), '-'), a.end());
    }
}

static void reverse_strings(Strings& aligned) {
    int size = aligned.size();
    for (int i = 0; i < size; i++) {
        std::string& a = aligned[i];
        std::reverse(a.begin(), a.end());
    }
}

static void process_seqs(Strings& seqs,
                         int mismatch_check,
                         int gap_check,
                         int aligned_check) {
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

static void append_seqs(Strings& aligned, const Strings& part) {
    int size = aligned.size();
    ASSERT_EQ(part.size(), size);
    for (int i = 0; i < size; i++) {
        aligned[i] += part[i];
    }
}

static int score_of(const Strings& rows) {
    Alignment aln((rows));
    int score = 0;
    int length = rows.front().length();
    for (int j = 0; j < length; j++) {
        if (is_equal(aln, j)) {
            score += 1;
        }
    }
    return score;
}

static void fix_bad_regions(Strings& aligned,
                            int mismatch_check,
                            int gap_check,
                            int aligned_check,
                            int min_length,
                            double min_identity) {
    Alignment aln((aligned));
    int length = aligned.front().length();
    std::vector<bool> good_col((length));
    for (int j = 0; j < length; j++) {
        good_col[j] = is_equal(aln, j);
    }
    int wf = FindLowSimilar::get_weight_factor(min_identity);
    Regions regions = FindLowSimilar::make_regions(good_col, wf);
    FindLowSimilar::reduce_regions(regions, min_length);
    Strings new_aligned((aln.size));
    BOOST_FOREACH (const Region& region, regions) {
        if (region.good_) {
            append_region(new_aligned, aln, region);
        } else {
            Strings seqs((aln.size));
            append_region(seqs, aln, region);
            int before_score = score_of(seqs);
            filter_out_gaps(seqs);
            reverse_strings(seqs);
            process_seqs(seqs, mismatch_check,
                         gap_check, aligned_check);
            int after_score = score_of(seqs);
            if (after_score > before_score) {
                reverse_strings(seqs);
                append_seqs(new_aligned, seqs);
            } else {
                append_region(new_aligned, aln, region);
            }
        }
    }
    aligned.swap(new_aligned);
}

static void realing_end(Strings& aligned,
                        int mismatch_check,
                        int gap_check,
                        int aligned_check) {
    int size = aligned.size();
    int length = aligned.front().length();
    if (length < 2) {
        return;
    }
    int prefix_length = length - aligned_check;
    if (prefix_length < 1) {
        prefix_length = 1;
    }
    Strings tails((size));
    for (int i = 0; i < size; i++) {
        std::string& row = aligned[i];
        tails[i] = row.substr(prefix_length);
        row.resize(prefix_length);
    }
    filter_out_gaps(tails);
    reverse_strings(tails);
    process_seqs(tails, mismatch_check,
                 gap_check, aligned_check);
    reverse_strings(tails);
    for (int i = 0; i < size; i++) {
        aligned[i] += tails[i];
    }
}

void SimilarAligner::similar_aligner(Strings& seqs) const {
    TimeIncrementer ti(this);
    if (seqs.empty()) {
        return;
    }
    int mismatch_check = opt_value("mismatch-check").as<int>();
    int gap_check = opt_value("gap-check").as<int>();
    int aligned_check = opt_value("aligned-check").as<int>();
    int min_length = opt_value("min-length").as<int>();
    double min_id = opt_value("min-identity").as<double>();
    process_seqs(seqs, mismatch_check, gap_check,
                 aligned_check);
    fix_bad_regions(seqs, mismatch_check,
                    gap_check, aligned_check,
                    min_length, min_id);
    realing_end(seqs, mismatch_check,
                gap_check, aligned_check);
}

SimilarAligner::SimilarAligner() {
    add_gopt("mismatch-check",
             "Min number of equal columns after single mismatch",
             "MISMATCH_CHECK");
    add_gopt("gap-check",
             "Min number of equal columns after single gap",
             "GAP_CHECK");
    add_gopt("aligned-check", "Min equal aligned part",
             "ALIGNED_CHECK");
    add_gopt("min-length", "Min length of fragment",
             "MIN_LENGTH");
    add_gopt("min-identity", "Min identity of block",
             "MIN_IDENTITY");
}

void SimilarAligner::align_seqs_impl(Strings& seqs) const {
    similar_aligner(seqs);
}

const char* SimilarAligner::name_impl() const {
    return "Align blocks with high similarity";
}

}

