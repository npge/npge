/*
 * NPG-explorer, Nucleotide PanGenome explorer
 * Copyright (C) 2012-2016 Boris Nagaev
 *
 * See the LICENSE file for terms of use.
 */

#ifndef NPGE_GENERAL_ALIGNER_HPP
#define NPGE_GENERAL_ALIGNER_HPP

#include <vector>
#include <algorithm>
#include <utility>

#include "global.hpp"
#include "throw_assert.hpp"
#include "Exception.hpp"
#include "cast.hpp"

namespace npge {

// TODO: gap_open

const int BAD_VALUE = 1e6;
const int MATRICES_NUMBER = 2;

/** Find the end of good alignment using Needleman-Wunsch with gap frame */
template <typename Contents>
class GeneralAligner {
public:
    enum Track {
        MATCH = 0,
        ROW_INC = +1,
        COL_INC = -1,
        STOP = +2
    };

    /** Constructor */
    GeneralAligner():
        gap_range_(1), max_errors_(0), gap_penalty_(1), local_(false) {
    }

    /** Get contents */
    const Contents& contents() const {
        return contents_;
    }

    /** Set contents
    Contents must have fillowing methods:
     - int first_size
     - int second_size
     - int substitution(pos_in_first, pos_in_second)
    */
    void set_contents(const Contents& contents) {
        contents_ = contents;
    }

    /** Get gap range */
    int gap_range() const {
        return gap_range_;
    }

    /** Set gap range.
    \code gap frame = 2 * gap range + 1 \endcode
    Gap range is max distance from main diagonal of considered
    states of pair alignment. The more gap_range, the more time.
    */
    void set_gap_range(int gap_range) {
        gap_range_ = gap_range;
    }

    /** Get max errors */
    int max_errors() const {
        return max_errors_;
    }

    /** Set max errors.
    Mismatch or gap are considered as 1 error.
    Alignment stops, when max errors accumulated.
    max_errors = -1 means no limit on errors.
    */
    void set_max_errors(int max_errors) {
        max_errors_ = max_errors;
    }

    /** Get gap penalty */
    int gap_penalty() const {
        return gap_penalty_;
    }

    /** Set gap penalty */
    void set_gap_penalty(int gap_penalty) {
        gap_penalty_ = gap_penalty;
    }

    /** Return if the alignment is local */
    bool local() const {
        return local_;
    }

    /** Set if the alignment is local.
    If local alignment is used, then max_errors() must be -1.
    Default: false.
    */
    void set_local(bool local) {
        local_ = local;
    }

    /** Run alignment algorithm.
    \param first_last Last aligned position in first sequence (output)
    \param second_last Last aligned position in second sequence (output)
    */
    void align(int& first_last, int& second_last) const {
        adjust_matrix_size();
        limit_range();
        make_frame();
        ASSERT_TRUE(!local() || max_errors() == -1);
        int& r_row = first_last;
        int& r_col = second_last;
        r_row = r_col = -1;
        for (int row = 0; row <= max_row(); row++) {
            int start_col = min_col(row);
            int stop_col = max_col(row);
            int min_score_col = start_col;
            for (int col = start_col; col <= stop_col; col++) {
                ASSERT_TRUE(col >= 0 && col < side());
                ASSERT_TRUE(in(row, col));
                int match = at(row - 1, col - 1) +
                            substitution(row, col);
                int gap1 = at(row, col - 1) + gap_penalty();
                int gap2 = at(row - 1, col) + gap_penalty();
                int score = std::min(match, std::min(gap1, gap2));
                if (local()) {
                    score = std::min(score, 0);
                }
                at(row, col) = score;
                if (score < at(row, min_score_col)) {
                    min_score_col = col;
                }
                track(row, col) = (score == match) ? MATCH :
                                  (score == gap1) ? COL_INC :
                                  ROW_INC;
            }
            if (max_errors() != -1 &&
                    at(row, min_score_col) > max_errors()) {
                break;
            }
            r_row = row;
            r_col = min_score_col;
        }
        if (max_errors() == -1) {
            // col -> max_col in this row
            ASSERT_TRUE(in(max_row(), max_col(max_row())));
            ASSERT_EQ(r_row, max_row());
            r_col = max_col(max_row());
            // if stopped earlier because of gap range
            int last_row = contents().first_size() - 1;
            int last_col = contents().second_size() - 1;
            if (r_row == last_row) {
                while (r_col < last_col) {
                    r_col += 1;
                    track(r_row, r_col) = COL_INC;
                }
            } else if (r_col == last_col) {
                while (r_row < last_row) {
                    r_row += 1;
                    track(r_row, r_col) = ROW_INC;
                }
            } else {
                throw Exception("row and column are not last");
            }
        }
    }

    /** Finds minimum cell <= (row, col) */
    void find_opt(int& row, int& col) const {
        int row0 = row;
        int col0 = col;
        for (int i = 0; i <= row0; i++) {
            for (int j = 0; j <= col0; j++) {
                if (at(i, j) < at(row, col)) {
                    row = i;
                    col = j;
                }
            }
        }
    }

    /** Return minimum value of matrix */
    int opt_score() const {
        int min_row = rows() - 1;
        int min_col = cols() - 1;
        find_opt(min_row, min_col);
        return at(min_row, min_col);
    }

    void local_to_global() {
        ASSERT_EQ(max_errors(), -1);
        track_local(rows() - 1, cols() - 1);
    }

    // ignores gap range
    void track_local(int row, int col) const {
        if (row == -1 || col == -1) {
            return;
        }
        int min_row = row;
        int min_col = col;
        find_opt(min_row, min_col);
        if (at(min_row, min_col) == 0) {
            return;
        }
        // go right to col
        for (int j = min_col; j <= col; j++) {
            track(min_row, j) = COL_INC;
        }
        // go bottom to row
        for (int i = min_row; i <= row; i++) {
            track(i, col) = ROW_INC;
        }
        while (at(min_row, min_col) < 0) {
            go_prev(min_row, min_col);
            ASSERT_TRUE(in(min_row, min_col));
        }
        track_local(min_row, min_col);
    }

    /** Strip out bad alignment tail.
    Consider the following alignment:
    \verbatim
    TTCCGGTGCTGCGaggga
    TTCCGGTGCTGCGcctct
    \endverbatim
    After tail stripping, it would be changed to
    \verbatim
    TTCCGGTGCTGCG
    TTCCGGTGCTGCG
    \endverbatim
    */
    void cut_tail(int& first_last, int& second_last) const {
        int& r_row = first_last;
        int& r_col = second_last;
        while (true) {
            int prev_row = r_row;
            int prev_col = r_col;
            go_prev(prev_row, prev_col);
            if (in(prev_row, prev_col) &&
                    at(prev_row, prev_col) < at(r_row, r_col)) {
                r_row = prev_row;
                r_col = prev_col;
            } else {
                break;
            }
        }
    }

    /** Write alignment as list of pairs of indices.
    \param first_last Last aligned position in first sequence (input)
    \param second_last Last aligned position in second sequence (input)
    \param alignment Array of pairs of positions, -1 as gaps
    */
    void export_alignment(int first_last, int second_last,
                          PairAlignment& alignment) const {
        int row = first_last, col = second_last;
        while (row != -1 || col != -1) {
            int tr = track(row, col);
            if (tr == STOP) {
                tr = MATCH;
            }
            bool print_first = (tr == MATCH || tr == ROW_INC);
            bool print_second = (tr == MATCH || tr == COL_INC);
            int a_row = print_first ? row : -1;
            int a_col = print_second ? col : -1;
            alignment.push_back(std::make_pair(a_row, a_col));
            if (track(row, col) == STOP) {
                break;
            }
            go_prev(row, col);
            ASSERT_TRUE(in(row, col));
        }
        std::reverse(alignment.begin(), alignment.end());
    }

    int rows() const {
        return contents_.first_size();
    }

    int cols() const {
        return contents_.second_size();
    }

    int rows_1() const {
        return rows() + 1;
    }

    int cols_1() const {
        return cols() + 1;
    }

    int side() const {
        return std::min(std::min(rows(), cols()) + gap_range(),
                        std::max(rows(), cols()));
    }

    int max_row() const {
        return std::min(rows(), side()) - 1;
    }

    int min_col(int row) const {
        return std::max(0, row - gap_range());
    }

    int max_col(int row) const {
        return std::min(cols() - 1, std::min(side() - 1,
                                             row + gap_range()));
    }

    int& at(int row0, int col0) const {
        ASSERT_MSG(in(row0, col0),
                   (TO_S(row0) + " " + TO_S(col0)).c_str());
        int row = row0 + 1, col = col0 + 1;
        int index = row * cols_1() + col;
        return matrix_[index * MATRICES_NUMBER];
    }

    /** Matrix with back track of alignment.
    \see Track
    */
    int& track(int row0, int col0) const {
        ASSERT_MSG(in(row0, col0),
                   (TO_S(row0) + " " + TO_S(col0)).c_str());
        int row = row0 + 1, col = col0 + 1;
        int index = row * cols_1() + col;
        return matrix_[index * MATRICES_NUMBER + 1];
    }

    /** Go to previous cell using track() */
    void go_prev(int& row, int& col) const {
        int tr = track(row, col);
        if (tr == MATCH || tr == ROW_INC) {
            row -= 1;
        }
        if (tr == MATCH || tr == COL_INC) {
            col -= 1;
        }
    }

    /** Go prev while at < 0, mark end with STOP */
    void find_stop(int& min_row, int& min_col) const {
        while (track(min_row, min_col) != STOP) {
            if (at(min_row, min_col) >= 0 ||
                    min_row == 0 || min_col == 0) {
                track(min_row, min_col) = STOP;
                break;
            }
            go_prev(min_row, min_col);
        }
    }

    bool in(int row, int col) const {
        bool row_is_good = (-1 <= row && row < rows());
        bool col_is_good = (-1 <= col && col < cols());
        return row_is_good && col_is_good;
    }

    int substitution(int row, int col) const {
        ASSERT_TRUE(row >= 0 && row < rows());
        ASSERT_TRUE(col >= 0 && col < cols());
        return contents_.substitution(row, col);
    }

    void adjust_matrix_size() const {
        int size = rows_1() * cols_1() * MATRICES_NUMBER;
        matrix_.resize(size, BAD_VALUE);
    }

    void limit_range() const {
        for (int row = -1; row < rows(); row++) {
            for (int ori = -1; ori <= 1; ori += 2) {
                int col = row + ori * (gap_range() + 1);
                if (in(row, col)) {
                    at(row, col) = BAD_VALUE;
                }
            }
        }
    }

    void make_frame() const {
        at(-1, -1) = 0;
        track(-1, -1) = STOP;
        for (int row = 0; row < rows(); row++) {
            if (local()) {
                at(row, -1) = 0;
            } else {
                at(row, -1) = (row + 1) * gap_penalty();
            }
            track(row, -1) = ROW_INC;
        }
        for (int col = 0; col < cols(); col++) {
            if (local()) {
                at(-1, col) = 0;
            } else {
                at(-1, col) = (col + 1) * gap_penalty();
            }
            track(-1, col) = COL_INC;
        }
    }

private:
    mutable std::vector<int> matrix_;
    int gap_range_, max_errors_, gap_penalty_;
    bool local_;
    Contents contents_;
};

template<typename Contents>
class ContentsProxy {
    Contents source_;
    int f_begin_;
    int f_length_;
    int f_ori_; // 1 or -1
    int s_begin_;
    int s_length_;
    int s_ori_; // 1 or -1

public:
    typedef ContentsProxy<Contents> Proxy;

    int first_size() const {
        return f_length_;
    }

    int second_size() const {
        return s_length_;
    }

    int substitution(int pos_in_first, int pos_in_second) const {
        ASSERT_LT(pos_in_first, f_length_);
        ASSERT_LT(pos_in_second, s_length_);
        int src_1 = map_f_to_src(pos_in_first);
        int src_2 = map_s_to_src(pos_in_second);
        return source_.substitution(src_1, src_2);
    }

    ContentsProxy() {
    }

    ContentsProxy(const Contents& source):
        source_(source),
        f_begin_(0),
        f_length_(source.first_size()),
        f_ori_(1),
        s_begin_(0),
        s_length_(source.second_size()),
        s_ori_(1) {
    }

    int map_f_to_src(int pos_in_first) const {
        int r = f_begin_ + pos_in_first * f_ori_;
        if (r < 0) {
            r += source_.first_size();
        }
        if (r >= source_.first_size()) {
            r -= source_.first_size();
        }
        return r;
    }

    int map_s_to_src(int pos_in_second) const {
        int r = s_begin_ + pos_in_second * s_ori_;
        if (r < 0) {
            r += source_.second_size();
        }
        if (r >= source_.second_size()) {
            r -= source_.second_size();
        }
        return r;
    }

    Proxy slice(int f_begin, int f_length,
                int s_begin, int s_length) const {
        ASSERT_GTE(f_length, 0);
        ASSERT_GTE(s_length, 0);
        Proxy result(source_);
        result.f_begin_ = map_f_to_src(f_begin);
        result.f_length_ = f_length;
        result.f_ori_ = f_ori_;
        result.s_begin_ = map_s_to_src(s_begin);
        result.s_length_ = s_length;
        result.s_ori_ = s_ori_;
        return result;
    }

    void alignment_to_src(PairAlignment& src,
                          const PairAlignment& dst) const {
        BOOST_FOREACH (const AlignmentPair& c, dst) {
            AlignmentPair src_c((-1), -1);
            if (c.first != -1) {
                src_c.first = map_f_to_src(c.first);
            }
            if (c.second != -1) {
                src_c.second = map_s_to_src(c.second);
            }
            src.push_back(src_c);
        }
    }
};

template <typename Proxy>
struct LocalAlignment {
    GeneralAligner<Proxy> ga_;
    int f_size_;
    int s_size_;
    int f_begin_;
    int s_begin_;
    int f_last_;
    int s_last_;
    int score_;

    int f_l() const {
        return f_begin_;
    }

    int s_l() const {
        return s_begin_;
    }

    int f_r() const {
        ASSERT_LT(f_last_, f_size_);
        return f_size_ - f_last_ - 1;
    }

    int s_r() const {
        ASSERT_LT(s_last_, s_size_);
        return s_size_ - s_last_ - 1;
    }

    const Proxy& proxy() const {
        return ga_.contents();
    }

    Proxy slice_left() const {
        return proxy().slice(0, f_l(), 0, s_l());
    }

    Proxy slice_right() const {
        return proxy().slice(f_last_ + 1, f_r(),
                             s_last_ + 1, s_r());
    }

    Proxy slice_both() const {
        int l1 = f_r() + f_l();
        int l2 = s_r() + s_l();
        return proxy().slice(f_last_ + 1, l1, s_last_ + 1, l2);
    }

    void export_src_aln(PairAlignment& aln) const {
        PairAlignment dst;
        ga_.export_alignment(f_last_, s_last_, dst);
        ASSERT_FALSE(dst.empty());
        ASSERT_TRUE(dst.front().first == f_begin_ ||
                    dst.front().second == s_begin_);
        ASSERT_TRUE(dst.back().first == f_last_ ||
                    dst.back().second == s_last_);
        proxy().alignment_to_src(aln, dst);
    }

    void export_dummy_src_aln(PairAlignment& aln) const {
        for (int i = 0; i < f_size_; i++) {
            int src_i = proxy().map_f_to_src(i);
            aln.push_back(std::make_pair(src_i, -1));
        }
        for (int j = 0; j < s_size_; j++) {
            int src_j = proxy().map_s_to_src(j);
            aln.push_back(std::make_pair(-1, src_j));
        }
    }

    LocalAlignment(const Proxy& contents, int gap_penalty) {
        ga_.set_max_errors(-1); // unlimited errors
        ga_.set_local(true);
        f_size_ = contents.first_size();
        s_size_ = contents.second_size();
        ga_.set_gap_penalty(gap_penalty);
        ga_.set_gap_range(std::max(f_size_, s_size_));
        ga_.set_contents(contents);
    }

    void build() {
        ASSERT_GT(f_size_, 0);
        ASSERT_GT(s_size_, 0);
        int _, __;
        ga_.align(_, __);
        f_last_ = f_size_ - 1;
        s_last_ = s_size_ - 1;
        ga_.find_opt(f_last_, s_last_);
        score_ = ga_.at(f_last_, s_last_);
        f_begin_ = f_last_;
        s_begin_ = s_last_;
        ga_.find_stop(f_begin_, s_begin_);
        ASSERT_GTE(f_begin_, 0);
        ASSERT_LT(f_begin_, f_size_);
        ASSERT_GTE(s_begin_, 0);
        ASSERT_LT(s_begin_, s_size_);
    }
};

template <typename Proxy>
int find_aln(PairAlignment& result, const Proxy& c,
             int gap_penalty, bool allow_shift) {
    LocalAlignment<Proxy> la((c), gap_penalty);
    if (la.f_size_ == 0 || la.s_size_ == 0) {
        la.export_dummy_src_aln(result);
        return 0;
    }
    la.build();
    int score = la.score_;
    if (score >= 0) {
        la.export_dummy_src_aln(result);
        return 0;
    }
    if (allow_shift) {
        la.export_src_aln(result);
        Proxy both = la.slice_both();
        score += find_aln(result, both, gap_penalty, false);
    } else {
        Proxy left = la.slice_left();
        score += find_aln(result, left, gap_penalty, false);
        //
        la.export_src_aln(result);
        //
        Proxy right = la.slice_right();
        score += find_aln(result, right, gap_penalty, false);
    }
    return score;
}

}

#endif

