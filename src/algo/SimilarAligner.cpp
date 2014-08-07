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
typedef std::vector<Ints> IntsCollection;
typedef std::map<int, int> Seq2Pos;
typedef std::map<std::string, Seq2Pos> Found;
typedef FindLowSimilar::Region Region;
typedef std::vector<Region> Regions;

struct Alignment {
    const Strings& seqs;
    Strings aligned;
    Ints pos;
    const int size;

    Alignment(const Strings& s):
        seqs(s), size(s.size()) {
        aligned.resize(size);
        pos.resize(size);
    }
};

struct SimilarAlignerImpl {
    int mismatch_check_;
    int gap_check_;
    int aligned_check_;
    int min_length_;
    Decimal min_identity_;

    bool equal_length(const Alignment& aln) const {
        int length = aln.aligned.front().length();
        for (int i = 1; i < aln.size; i++) {
            if (aln.aligned[i].length() != length) {
                return false;
            }
        }
        return true;
    }

    bool is_stop(const Alignment& aln, int shift = 0) const {
        for (int i = 0; i < aln.size; i++) {
            if (aln.pos[i] + shift >= aln.seqs[i].size()) {
                return true;
            }
        }
        return false;
    }

    void append_cols(Alignment& aln, int cols = 1) const {
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

    void append_gaps(Alignment& aln) const {
        int max_l = 0;
        for (int i = 0; i < aln.size; i++) {
            int l = aln.aligned[i].length();
            max_l = std::max(max_l, l);
        }
        for (int i = 0; i < aln.size; i++) {
            aln.aligned[i].resize(max_l, '-');
        }
        ASSERT_TRUE(equal_length(aln));
    }

    void append_all(Alignment& aln) const {
        for (int i = 0; i < aln.size; i++) {
            int& p = aln.pos[i];
            std::string tail = aln.seqs[i].substr(p);
            aln.aligned[i] += tail;
            p += tail.size();
        }
        append_gaps(aln);
    }

    bool is_equal(const Ints& pos, const Alignment& aln,
                  int shift = 0, int cols = 1) const {
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

    bool is_equal(const Alignment& aln,
                  int shift = 0, int cols = 1) const {
        return is_equal(aln.pos, aln, shift, cols);
    }

    bool is_mismatch(const Alignment& aln) const {
        return !is_stop(aln, mismatch_check_) &&
               is_equal(aln, 1, mismatch_check_);
    }

    bool try_mismatch(Alignment& aln) const {
        if (is_mismatch(aln)) {
            append_cols(aln, mismatch_check_ + 1);
            return true;
        } else {
            return false;
        }
    }

    bool append_chars(Alignment& aln, int i,
                      int cols = 1) const {
        int& p = aln.pos[i];
        for (int j = 0; j < cols; j++) {
            aln.aligned[i] += aln.seqs[i][p];
            p += 1;
        }
    }

    bool make_gap_shift(Ints& equal_pos, char c,
                        const Alignment& aln) const {
        equal_pos.resize(aln.size);
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
        int shift = 0;
        return is_equal(equal_pos, aln, shift, gap_check_);
    }

    void apply_gap(Alignment& aln, const Ints& equal_pos,
                   int gap_check) const {
        for (int i = 0; i < aln.size; i++) {
            if (equal_pos[i] == aln.pos[i] + 1) {
                append_chars(aln, i);
            }
        }
        append_gaps(aln);
        append_cols(aln, gap_check);
    }

    void find_all_gaps(IntsCollection& variants,
                       const Alignment& aln) const {
        std::set<char> chars;
        for (int i = 0; i < aln.size; i++) {
            int p = aln.pos[i];
            chars.insert(aln.seqs[i][p]);
        }
        BOOST_FOREACH (char c, chars) {
            Ints equal_pos;
            if (make_gap_shift(equal_pos, c, aln)) {
                variants.push_back(equal_pos);
            }
        }
    }

    void find_best_gap(IntsCollection& variants,
                       Alignment& aln) const {
        // find the best of them by increasing gap_check
        for (int gap_check = gap_check_ + 1;; gap_check += 1) {
            IntsCollection next_variants;
            BOOST_FOREACH (const Ints& equal_pos, variants) {
                int shift = 0;
                if (is_equal(equal_pos, aln,
                            shift, gap_check)) {
                    next_variants.push_back(equal_pos);
                }
            }
            if (next_variants.empty()) {
                // can not find the best variant
                // use one of variants for previous gap_check
                apply_gap(aln, variants.front(), gap_check - 1);
                return;
            } else if (next_variants.size() == 1) {
                // the best variant was found
                apply_gap(aln, next_variants.front(), gap_check);
                return;
            } else {
                // go on
                variants.swap(next_variants);
            }
        }
    }

    bool try_gap(Alignment& aln) const {
        if (is_stop(aln, gap_check_)) {
            return false;
        }
        IntsCollection variants;
        find_all_gaps(variants, aln);
        if (variants.empty()) {
            return false;
        }
        if (variants.size() == 1) {
            apply_gap(aln, variants.front(), gap_check_);
            return true;
        }
        // several variants
        find_best_gap(variants, aln);
        return true;
    }

    int min_tail(const Alignment& aln) const {
        int mt = aln.seqs.front().size() - aln.pos.front();
        for (int i = 1; i < aln.size; i++) {
            int t = aln.seqs[i].size() - aln.pos[i];
            mt = std::min(mt, t);
        }
        return mt;
    }

    std::string find_best_word(const Alignment& aln,
                               Found& ff, int shift) const {
        typedef std::set<std::string> StringSet;
        StringSet words;
        std::string best_word;
        for (int i = 0; i < aln.size; i++) {
            int p = aln.pos[i] + shift;
            const std::string& seq = aln.seqs[i];
            std::string word = seq.substr(p, aligned_check_);
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

    void append_aligned(Alignment& aln,
                        const Seq2Pos& s2p) const {
        Strings tmp_seqs((aln.size));
        for (int i = 0; i < aln.size; i++) {
            int p = aln.pos[i];
            int length = s2p.find(i)->second;
            std::string& tmp_seq = tmp_seqs[i];
            tmp_seq = aln.seqs[i].substr(p, length);
            std::reverse(tmp_seq.begin(), tmp_seq.end());
        }
        process_seqs(tmp_seqs);
        for (int i = 0; i < aln.size; i++) {
            std::string& tmp_row = tmp_seqs[i];
            std::reverse(tmp_row.begin(), tmp_row.end());
            aln.aligned[i] += tmp_row;
            int& p = aln.pos[i];
            int length = s2p.find(i)->second;
            p += length;
        }
    }

    bool try_aligned(Alignment& aln) const {
        Found ff;
        int max_shift = min_tail(aln) - aligned_check_;
        for (int shift = 0; shift < max_shift; shift++) {
            std::string best_word = find_best_word(aln,
                    ff, shift);
            if (!best_word.empty()) {
                append_aligned(aln, ff[best_word]);
                append_cols(aln, aligned_check_);
                return true;
            }
        }
        return false;
    }

    bool pos_less(const Ints& a, const Ints& b,
                  int shift = 0) const {
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
    void append_end(Alignment& aln) const {
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

    void process_cols(Alignment& aln) const {
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

    void append_region(Strings& aligned, const Alignment& aln,
                       const Region& r) const {
        for (int i = 0; i < aln.size; i++) {
            const std::string& seq = aln.seqs[i];
            aligned[i] += seq.substr(r.start_, r.length());
        }
    }

    void filter_out_gaps(Strings& aligned) const {
        int size = aligned.size();
        for (int i = 0; i < size; i++) {
            std::string& a = aligned[i];
            a.erase(std::remove(a.begin(), a.end(), '-'),
                    a.end());
        }
    }

    void reverse_strings(Strings& aligned) const {
        int size = aligned.size();
        for (int i = 0; i < size; i++) {
            std::string& a = aligned[i];
            std::reverse(a.begin(), a.end());
        }
    }

    void process_seqs(Strings& seqs) const {
        Alignment aln((seqs));
        process_cols(aln);
        for (int i = 0; i < aln.size; i++) {
            ASSERT_GTE(aln.aligned[i].size(),
                       aln.seqs[i].size());
        }
        ASSERT_TRUE(equal_length(aln));
        seqs.swap(aln.aligned);
    }

    void append_seqs(Strings& aligned,
                     const Strings& part) const {
        int size = aligned.size();
        ASSERT_EQ(part.size(), size);
        for (int i = 0; i < size; i++) {
            aligned[i] += part[i];
        }
    }

    int score_of(const Strings& rows) const {
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

    void fix_bad_regions(Strings& aligned) const {
        Alignment aln((aligned));
        int length = aligned.front().length();
        std::vector<bool> good_col((length));
        for (int j = 0; j < length; j++) {
            good_col[j] = is_equal(aln, j);
        }
        int wf = FindLowSimilar::get_weight_factor(min_identity_);
        Regions regions = FindLowSimilar::make_regions(good_col, wf);
        FindLowSimilar::reduce_regions(regions, min_length_);
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
                process_seqs(seqs);
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

    void realing_end(Strings& aligned) const {
        int size = aligned.size();
        int length = aligned.front().length();
        if (length < 2) {
            return;
        }
        int prefix_length = length - aligned_check_;
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
        process_seqs(tails);
        reverse_strings(tails);
        for (int i = 0; i < size; i++) {
            aligned[i] += tails[i];
        }
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
    im.min_length_ = opt_value("min-length").as<int>();
    im.min_identity_ = opt_value("min-identity").as<Decimal>();
    im.process_seqs(seqs);
    im.fix_bad_regions(seqs);
    im.realing_end(seqs);
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

std::string SimilarAligner::aligner_type() const {
    return "similar";
}

const char* SimilarAligner::name_impl() const {
    return "Align blocks with high similarity";
}

}

