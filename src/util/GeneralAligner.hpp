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

/** Find the end of good alignment using Needleman-Wunsch with gap frame */
template <typename Contents>
class GeneralAligner {
public:
    /** Pair alignment, array of pairs of positions, -1 as gaps */
    typedef std::vector<std::pair<int, int> > PairAlignment;

    /** Constructor */
    GeneralAligner():
        gap_range_(1), max_errors_(0), gap_penalty_(1)
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
    */
    void set_max_errors(int max_errors) {
        max_errors_ = max_errors;
        gap_range_ = std::min(gap_range_, max_errors_ / gap_penalty_);
    }

    /** Get gap penalty */
    int gap_penalty() const {
        return gap_penalty_;
    }

    /** Set gap penalty */
    void set_gap_penalty(int gap_penalty) {
        gap_penalty_ = gap_penalty;
    }

    /** Run alignment algorithm.
    \param first_last Last aligned position in first sequence (output)
    \param second_last Last aligned position in second sequence (output)
    */
    void align(int& first_last, int& second_last) const {
        adjust_matrix_size();
        int& r_row = first_last;
        int& r_col = second_last;
        r_row = r_col = -1;
        if (in(0, 0)) {
            at(0, 0) = 0;
        }
        for (int row = 0; row <= max_row(); row++) {
            int start_col = min_col(row);
            int min_score_col = start_col;
            for (int col = start_col; col <= max_col(row); col++) {
                BOOST_ASSERT(col >= 0 && col < side());
                BOOST_ASSERT(in(row, col));
                int score = substitution(row, col);
                if (in(row - 1, col - 1)) {
                    score += at(row - 1, col - 1);
                }
                if (in(row - 1, col)) {
                    int alt_score = at(row - 1, col) + gap_penalty();
                    score = std::min(score, alt_score);
                }
                if (in(row, col - 1)) {
                    int alt_score = at(row, col - 1) + gap_penalty();
                    score = std::min(score, alt_score);
                }
                at(row, col) = score;
                if (score < at(row, min_score_col)) {
                    min_score_col = col;
                }
            }
            if (at(row, min_score_col) > max_errors()) {
                break;
            }
            r_row = row;
            r_col = min_score_col;
        }
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
            if (in(r_row - 1, r_col) &&
                    at(r_row - 1, r_col) < at(r_row, r_col)) {
                r_row -= 1;
            } else if (in(r_row, r_col - 1) &&
                       at(r_row, r_col - 1) < at(r_row, r_col)) {
                r_col -= 1;
            } else if (in(r_row - 1, r_col - 1) &&
                       at(r_row - 1, r_col - 1) < at(r_row, r_col)) {
                r_row -= 1;
                r_col -= 1;
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
        while (row >= 0 && col >= 0) {
            bool print_first = true;
            bool print_second = true;
            if (in(row - 1, col) && at(row - 1, col) < at(row, col)) {
                print_second = false;
            } else if (in(row, col - 1) &&
                       at(row, col - 1) < at(row, col)) {
                print_first = false;
            }
            int a_row = print_first ? row : -1;
            int a_col = print_second ? col : -1;
            alignment.push_back(std::make_pair(a_row, a_col));
            if (print_first) {
                row -= 1;
            }
            if (print_second) {
                col -= 1;
            }
        }
        std::reverse(alignment.begin(), alignment.end());
    }

    int rows() const {
        return contents_.first_size();
    }

    int cols() const {
        return contents_.second_size();
    }

    int side() const {
        return std::min(std::min(rows(), cols()) + gap_range(),
                        std::max(rows(), cols()));
    }

    int row_size() const {
        return 1 + 2 * gap_range();
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

    int& at(int row, int col) const {
        return matrix_[row * row_size() + gap_range() + col - row];
    }

    bool in(int row, int col) const {
        return row >= 0 && row < side() &&
               col >= min_col(row) && col <= max_col(row);
    }

    int substitution(int row, int col) const {
        return contents_.substitution(row, col);
    }

    void adjust_matrix_size() const {
        matrix_.resize(side() * row_size());
    }

private:
    mutable std::vector<int> matrix_;
    int gap_range_, max_errors_, gap_penalty_;
    Contents contents_;
};

}

#endif

