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
#include "throw_assert.hpp"
#include "global.hpp"

namespace npge {

static const std::string empty_string_;

class Slice {
public:
    Slice():
        s_(&empty_string_),
        start_(0), length_(0), ori_(1) {
    }

    Slice(const Slice& other):
        s_(other.s_),
        start_(other.start_),
        length_(other.length_),
        ori_(other.ori_) {
    }

    Slice(const std::string& s):
        s_(&s), start_(0), length_(s.length()), ori_(1) {
    }

    int length() const {
        return length_;
    }

    bool empty() const {
        return length() == 0;
    }

    char at(int i) const {
        return (*s_)[source_index(i)];
    }

    void inverse() {
        start_ = source_index(length() - 1);
        ori_ *= -1;
    }

    void cut(int start, int length) {
        start_ = source_index(start);
        length_ = length;
    }

    void cut(int start) {
        cut(start, length() - start);
    }

    std::string to_s() const {
        std::string result;
        int l = length();
        result.reserve(l);
        for (int i = 0; i < l; i++) {
            result.push_back(at(i));
        }
        return result;
    }

    bool operator==(const Slice& other) const {
        if (length() != other.length()) {
            return false;
        }
        for (int i = 0; i < length(); i++) {
            if (at(i) != other.at(i)) {
                return false;
            }
        }
        return true;
    }

    bool operator!=(const Slice& other) const {
        return !(*this == other);
    }

    bool operator<(const Slice& other) const {
        int L = std::min(length(), other.length());
        for (int i = 0; i < L; i++) {
            if (at(i) < other.at(i)) {
                return true;
            } else if (at(i) > other.at(i)) {
                return false;
            }
        }
        if (length() < other.length()) {
            return true;
        } else {
            return false;
        }
    }

private:
    const std::string* s_;
    int start_;
    int length_;
    int ori_;

    int source_index(int i) const {
        return start_ + i * ori_;
    }
};

class Slices : public std::vector<Slice> {
public:
    Slices() {
    }

    Slices(const Strings& strings) {
        reserve(strings.size());
        BOOST_FOREACH (const std::string& s, strings) {
            push_back(Slice(s));
        }
    }

    void inverse() {
        BOOST_FOREACH (Slice& s, *this) {
            s.inverse();
        }
    }

    void cut(int start) {
        BOOST_FOREACH (Slice& s, *this) {
            s.cut(start);
        }
    }

    bool empty() const {
        BOOST_FOREACH (const Slice& s, *this) {
            if (!s.empty()) {
                return false;
            }
        }
        return true;
    }
};

typedef std::map<const Slice*, int> Seq2Pos;
typedef std::map<Slice, Seq2Pos> Found;

struct SimilarAlignerImpl {
    int mismatch_check_;
    int gap_check_;
    int aligned_check_;

    bool equal_length(const Strings& aligned) const {
        ASSERT_GT(aligned.size(), 0);
        int length = aligned.front().length();
        BOOST_FOREACH (const std::string& row, aligned) {
            if (row.length() != length) {
                return false;
            }
        }
        return true;
    }

    void append_gaps(Strings& aligned) const {
        int max_l = 0;
        BOOST_FOREACH (const std::string& row, aligned) {
            int l = row.length();
            max_l = std::max(max_l, l);
        }
        BOOST_FOREACH (std::string& row, aligned) {
            row.resize(max_l, '-');
        }
        ASSERT_TRUE(equal_length(aligned));
    }

    void take_short_seqs(Slices& long_seqs, Slices& short_seqs,
                         const Slices& seqs0) const {
        ASSERT_TRUE(long_seqs.empty());
        ASSERT_TRUE(short_seqs.empty());
        ASSERT_FALSE(seqs0.empty());
        int l = min_length(seqs0);
        BOOST_FOREACH (const Slice& seq, seqs0) {
            if (seq.length() > l) {
                long_seqs.push_back(seq);
            } else {
                short_seqs.push_back(seq);
            }
        }
    }

    void back_short_seqs(Strings& aligned0,
                         Strings& long_aligned,
                         Strings& short_aligned,
                         const Slices& seqs0) const {
        ASSERT_TRUE(aligned0.empty());
        ASSERT_FALSE(seqs0.empty());
        ASSERT_LTE(long_aligned.size(), seqs0.size());
        ASSERT_LTE(short_aligned.size(), seqs0.size());
        int l = min_length(seqs0);
        int long_index = 0;
        int short_index = 0;
        aligned0.resize(seqs0.size());
        for (int i = 0; i < seqs0.size(); i++) {
            const Slice& seq = seqs0[i];
            if (seq.length() > l) {
                aligned0[i].swap(long_aligned[long_index]);
                long_index += 1;
            } else {
                aligned0[i].swap(short_aligned[short_index]);
                short_index += 1;
            }
        }
        ASSERT_EQ(long_index, long_aligned.size());
        ASSERT_EQ(short_index, short_aligned.size());
        append_gaps(aligned0);
    }

    bool is_stop(const Slices& seqs, int min_size = 1) const {
        BOOST_FOREACH (const Slice& seq, seqs) {
            if (seq.length() < min_size) {
                return true;
            }
        }
        return false;
    }

    bool is_equal(const Slices& seqs, int cols = 1) const {
        if (is_stop(seqs, cols)) {
            return false;
        }
        for (int j = 0; j < cols; j += 1) {
            char c = seqs.front().at(j);
            BOOST_FOREACH (const Slice& seq, seqs) {
                if (seq.at(j) != c) {
                    return false;
                }
            }
        }
        return true;
    }

    void move_cols(Strings& aligned, Slices& seqs,
                   int cols = 1) const {
        ASSERT_EQ(aligned.size(), seqs.size());
        for (int col = 0; col < cols; col++) {
            int size = seqs.size();
            for (int i = 0; i < size; i++) {
                Slice& seq = seqs[i];
                std::string& row = aligned[i];
                ASSERT_GT(seq.length(), 0);
                row.push_back(seq.at(0));
                seq.cut(1);
            }
        }
    }

    bool is_mismatch(const Slices& seqs) const {
        if (is_stop(seqs, mismatch_check_ + 1)) {
            return false;
        }
        Slices copy = seqs;
        BOOST_FOREACH (Slice& seq, copy) {
            seq.cut(1);
        }
        return is_equal(copy, mismatch_check_);
    }

    typedef std::set<char> Chars;
    typedef std::map<char, Slices> Variants;

    void make_gap_variant(Slices& variant, char gap_char,
                          const Slices& seqs) const {
        BOOST_FOREACH (Slice seq, seqs) {
            ASSERT_GT(seq.length(), 0);
            if (seq.at(0) == gap_char) {
                seq.cut(1);
            }
            variant.push_back(seq);
        }
    }

    void select_best_variant(Variants& variants) const {
        while (variants.size() > 1) {
            Chars bad_chars;
            BOOST_FOREACH (Variants::value_type& v, variants) {
                char gap_char = v.first;
                Slices& seqs = v.second;
                if (!is_equal(seqs)) {
                    bad_chars.insert(gap_char);
                }
            }
            BOOST_FOREACH (char bad_char, bad_chars) {
                if (variants.size() > 1) {
                    variants.erase(bad_char);
                }
            }
            BOOST_FOREACH (Variants::value_type& v, variants) {
                Slices& seqs = v.second;
                seqs.cut(1);
            }
        }
    }

    char best_gap_char(const Slices& seqs) const {
        Chars chars;
        BOOST_FOREACH (const Slice& seq, seqs) {
            chars.insert(seq.at(0));
        }
        Variants variants;
        BOOST_FOREACH (char gap_char, chars) {
            Slices& variant = variants[gap_char];
            make_gap_variant(variant, gap_char, seqs);
        }
        select_best_variant(variants);
        ASSERT_EQ(variants.size(), 1);
        char gap_char = variants.begin()->first;
        Slices best_variant;
        make_gap_variant(best_variant, gap_char, seqs);
        if (!is_equal(best_variant, gap_check_)) {
            return 0;
        }
        return gap_char;
    }

    bool try_gap(Strings& aligned, Slices& seqs) const {
        if (is_stop(seqs, gap_check_)) {
            return false;
        }
        char gap_char = best_gap_char(seqs);
        if (gap_char == 0) {
            return false;
        }
        int size = seqs.size();
        for (int i = 0; i < size; i++) {
            Slice& seq = seqs[i];
            std::string& row = aligned[i];
            if (seq.at(0) == gap_char) {
                row.push_back(gap_char);
                seq.cut(1);
            } else {
                row.push_back('-');
            }
        }
        move_cols(aligned, seqs, gap_check_);
        return true;
    }

    void move_good_alignment(Strings& aligned,
                             Slices& seqs) const {
        ASSERT_EQ(aligned.size(), seqs.size());
        while (true) {
            if (is_equal(seqs)) {
                move_cols(aligned, seqs);
            } else if (is_mismatch(seqs)) {
                move_cols(aligned, seqs, mismatch_check_ + 1);
            } else if (try_gap(aligned, seqs)) {
                // ok
            } else {
                return;
            }
            ASSERT_TRUE(equal_length(aligned));
        }
    }

    void move_perfect_alignment(Strings& aligned,
                                Slices& seqs) const {
        ASSERT_EQ(aligned.size(), seqs.size());
        while (is_equal(seqs)) {
            move_cols(aligned, seqs);
        }
    }

    int min_length(const Slices& seqs) const {
        int result = seqs.front().length();
        BOOST_FOREACH (const Slice& seq, seqs) {
            result = std::min(result, seq.length());
        }
        return result;
    }

    int max_length(const Slices& seqs) const {
        int result = seqs.front().length();
        BOOST_FOREACH (const Slice& seq, seqs) {
            result = std::max(result, seq.length());
        }
        return result;
    }

    Slice find_best_word(const Slices& seqs,
                         Found& ff, int shift) const {
        typedef std::set<Slice> SliceSet;
        SliceSet words;
        Slice best_word;
        BOOST_FOREACH (const Slice& seq, seqs) {
            Slice word = seq;
            word.cut(shift, aligned_check_);
            words.insert(word);
            Seq2Pos& s2p = ff[word];
            if (s2p.find(&seq) == s2p.end()) {
                s2p[&seq] = shift;
            }
            if (s2p.size() == seqs.size()) {
                best_word = word;
            }
        }
        if (words.size() == 1) {
            // same word with shift
            best_word = *words.begin();
            BOOST_FOREACH (const Slice& seq, seqs) {
                ff[best_word][&seq] = shift;
            }
        }
        return best_word;
    }

    void find_good_alignment(Slices& bad, Slices& good,
                             const Slices& seqs) const {
        bad.clear();
        good.clear();
        Found ff;
        int max_shift = min_length(seqs) - aligned_check_;
        for (int shift = 0; shift < max_shift; shift++) {
            Slice best_word = find_best_word(seqs, ff, shift);
            if (!best_word.empty()) {
                Seq2Pos& s2p = ff[best_word];
                BOOST_FOREACH (const Slice& seq, seqs) {
                    int pos = s2p[&seq];
                    bad.push_back(seq);
                    bad.back().cut(0, pos);
                    good.push_back(seq);
                    good.back().cut(pos);
                }
                return;
            }
        }
        // whole sequences are bad
        good.resize(seqs.size());
        bad = seqs;
    }

    void reverse_strings(Strings& aligned) const {
        BOOST_FOREACH (std::string& row, aligned) {
            std::reverse(row.begin(), row.end());
        }
    }

    void push_back_strings(Strings& a,
                           const Strings& b) const {
        ASSERT_EQ(a.size(), b.size());
        for (int i = 0; i < a.size(); i++) {
            a[i] += b[i];
        }
    }

    bool is_empty(const Slices& seqs) const {
        BOOST_FOREACH (const Slice& seq, seqs) {
            if (!seq.empty()) {
                return false;
            }
        }
        return true;
    }

    void process_slices(Strings& aligned, Slices& seqs) const {
        if (seqs.empty()) {
            return;
        }
        ASSERT_EQ(aligned.size(), seqs.size());
        BOOST_FOREACH (const std::string& row, aligned) {
            ASSERT_EQ(row.length(), 0);
        }
        int max_l = max_length(seqs);
        // perfect alignment from right end
        Strings perfect_right;
        perfect_right.resize(seqs.size());
        seqs.inverse();
        move_perfect_alignment(perfect_right, seqs);
        seqs.inverse();
        reverse_strings(perfect_right);
        while (!is_empty(seqs)) {
            move_good_alignment(aligned, seqs);
            Slices bad, good;
            find_good_alignment(bad, good, seqs);
            if (!is_empty(bad)) {
                bad.inverse();
                Strings bad_right_al;
                bad_right_al.resize(bad.size());
                move_good_alignment(bad_right_al, bad);
                reverse_strings(bad_right_al);
                bad.inverse();
                Slices bad1long, bad1short;
                take_short_seqs(bad1long, bad1short, bad);
                Strings bad_al;
                if (!is_empty(bad1long)) {
                    Strings bad1long_al, bad1short_al;
                    bad1long_al.resize(bad1long.size());
                    bad1short_al.resize(bad1short.size());
                    process_slices(bad1long_al, bad1long);
                    process_slices(bad1short_al, bad1short);
                    back_short_seqs(bad_al, bad1long_al,
                                    bad1short_al, bad);
                } else {
                    // all are short
                    BOOST_FOREACH (const Slice& seq, bad) {
                        bad_al.push_back(seq.to_s());
                    }
                }
                push_back_strings(aligned, bad_al);
                push_back_strings(aligned, bad_right_al);
            }
            seqs = good;
        }
        push_back_strings(aligned, perfect_right);
        BOOST_FOREACH (const std::string& row, aligned) {
            ASSERT_GTE(row.length(), max_l);
        }
    }

    void process_seqs(Strings& seqs) const {
        Slices slices(seqs);
        Strings aligned;
        aligned.resize(seqs.size());
        process_slices(aligned, slices);
        seqs.swap(aligned);
    }
};

void SimilarAligner::similar_aligner(Strings& seqs) const {
    TimeIncrementer ti(this);
    if (seqs.empty()) {
        return;
    }
    SimilarAlignerImpl im;
    im.mismatch_check_ = opt_value("mismatch-check").as<int>();
    im.gap_check_ = opt_value("gap-check").as<int>();
    im.aligned_check_ = opt_value("aligned-check").as<int>();
    im.process_seqs(seqs);
}

SimilarAligner::SimilarAligner() {
    add_gopt("mismatch-check",
             "Min number of equal columns after "
             "single mismatch", "MISMATCH_CHECK");
    add_gopt("gap-check",
             "Min number of equal columns after single gap",
             "GAP_CHECK");
    add_gopt("aligned-check", "Min equal aligned part",
             "ALIGNED_CHECK");
}

void SimilarAligner::align_seqs_impl(Strings& seqs) const {
    similar_aligner(seqs);
}

std::string SimilarAligner::aligner_type() const {
    return "similar";
}

const char* SimilarAligner::name_impl() const {
    return "Align blocks with high similarity";
}

}

