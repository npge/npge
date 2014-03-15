/*
 * bloomrepeats, Find genomic repeats, using Bloom filter based prefiltration
 * Copyright (C) 2014 Boris Nagaev
 *
 * See the LICENSE file for terms of use.
 */

#ifndef BR_GENERAL_ALIGNER_HPP
#define BR_GENERAL_ALIGNER_HPP

#include <vector>
#include <algorithm>
#include <utility>

#include "global.hpp"
#include "throw_assert.hpp"

namespace bloomrepeats {

// TODO: gap_open

const int BAD_VALUE = 1e6;
const int MATRICES_NUMBER = 2;

/** Find the end of good alignment using Needleman-Wunsch with gap frame */
template <typename Contents>
class GeneralAligner {
public:
    /** Pair alignment, array of pairs of positions, -1 as gaps */
    typedef std::vector<std::pair<int, int> > PairAlignment;

    /** Constructor */
    GeneralAligner():
        gap_range_(1), max_errors_(0), gap_penalty_(1), local_(false)
    { }

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
        BOOST_ASSERT(!local() || max_errors() == -1);
        int& r_row = first_last;
        int& r_col = second_last;
        r_row = r_col = -1;
        for (int row = 0; row <= max_row(); row++) {
            int start_col = min_col(row);
            int stop_col = max_col(row);
            int min_score_col = start_col;
            for (int col = start_col; col <= stop_col; col++) {
                BOOST_ASSERT(col >= 0 && col < side());
                BOOST_ASSERT(in(row, col));
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
                track(row, col) = (score == match) ? 0 :
                                  (score == gap1) ? (-1) :
                                  (+1);
            }
            if (max_errors() != -1 &&
                    at(row, min_score_col) > max_errors()) {
                break;
            }
            r_row = row;
            r_col = min_score_col;
        }
        if (local()) {
            BOOST_ASSERT(max_errors() == -1);
            track_local(rows() - 1, cols() - 1);
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
            track(min_row, j) = -1;
        }
        // go bottom to row
        for (int i = min_row; i <= row; i++) {
            track(i, col) = 1;
        }
        while (at(min_row, min_col) < 0) {
            go_prev(min_row, min_col);
            BOOST_ASSERT(in(min_row, min_col));
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
            bool print_first = (tr == 0 || tr == 1);
            bool print_second = (tr == 0 || tr == -1);
            int a_row = print_first ? row : -1;
            int a_col = print_second ? col : -1;
            alignment.push_back(std::make_pair(a_row, a_col));
            go_prev(row, col);
            BOOST_ASSERT(in(row, col));
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
        BOOST_ASSERT(in(row0, col0));
        int row = row0 + 1, col = col0 + 1;
        int index = row * cols_1() + col;
        return matrix_[index * MATRICES_NUMBER];
    }

    /** Matrix with back track of alignment.
    0 means match,
    +1 means row increment,
    -1 means column increment.
    */
    int& track(int row0, int col0) const {
        BOOST_ASSERT(in(row0, col0));
        int row = row0 + 1, col = col0 + 1;
        int index = row * cols_1() + col;
        return matrix_[index * MATRICES_NUMBER + 1];
    }

    /** Go to previous cell using track() */
    void go_prev(int& row, int& col) const {
        int tr = track(row, col);
        if (tr == 0 || tr == 1) {
            row -= 1;
        }
        if (tr == 0 || tr == -1) {
            col -= 1;
        }
    }

    bool in(int row, int col) const {
        bool row_is_good = (-1 <= row && row < rows());
        bool col_is_good = (-1 <= col && col < cols());
        return row_is_good && col_is_good;
    }

    int substitution(int row, int col) const {
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
        for (int row = 0; row < rows(); row++) {
            if (local()) {
                at(row, -1) = 0;
            } else {
                at(row, -1) = (row + 1) * gap_penalty();
            }
            track(row, -1) = 1; // row increment
        }
        for (int col = 0; col < cols(); col++) {
            if (local()) {
                at(-1, col) = 0;
            } else {
                at(-1, col) = (col + 1) * gap_penalty();
            }
            track(-1, col) = -1; // col increment
        }
    }

private:
    mutable std::vector<int> matrix_;
    int gap_range_, max_errors_, gap_penalty_;
    bool local_;
    Contents contents_;
};

}

#endif

