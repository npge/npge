/*
 * bloomrepeats, Find genomic repeats, using Bloom filter based prefiltration
 * Copyright (C) 2012 Boris Nagaev
 *
 * See the LICENSE file for terms of use.
 */

#ifndef BR_PAIR_ALIGNER_HPP
#define BR_PAIR_ALIGNER_HPP

#include <vector>
#include <string>

#include "global.hpp"

namespace bloomrepeats {

/** Find the end of good alignment using Needleman-Wunsch with gap frame */
class PairAligner {
public:
    /** Constructor.
    \param max_errors Max number of errors in pair alignment.
        Alignment stops, when max errors accumulated.
    \param gap_range Max distance from main diagonal of considered
        states of pair alignment. The more gap_range, the more time.
    */
    PairAligner(int max_errors = 10, int gap_range = 5);

    /** Set first sequence */
    void set_first(const char* start, int size);

    /** Set second sequence */
    void set_second(const char* start, int size);

    /** Get gap range */
    int gap_range() const {
        return gap_range_;
    }

    /** Set gap range.
    \code gap frame = 2 * gap range + 1 \endcode
    Gap range is max distance from main diagonal of considered
    states of pair alignment. The more gap_range, the more time.
    */
    void set_gap_range(int gap_range);

    /** Get max errors */
    int max_errors() const {
        return max_errors_;
    }

    /** Set max errors.
    Mismatch or gap are considered as 1 error.
    Alignment stops, when max errors accumulated.
    */
    void set_max_errors(int max_errors);

    /** Get whether bad alignment tail would be stripped out */
    bool no_tail() const {
        return no_tail_;
    }

    /** Set whether bad alignment tail would be stripped out.
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

    Defaults to \c true.
    */
    void set_no_tail(bool no_tail) {
        no_tail_ = no_tail;
    }

    /** Run alignment algorithm.
    \param first_last Last aligned position in first sequence (output)
    \param second_last Last aligned position in second sequence (output)
    \param first_str Alignment row (like "a-tgc") (optional output)
    \param second_str Alignment row (like "a-tgc") (optional output)
    \param alignment Array of pairs of positions, -1 as gaps (optional output)
    \param gap Char used as gap for first_str and second_str
    */
    void align(int& first_last, int& second_last,
               std::string* first_str = 0, std::string* second_str = 0,
               std::vector<std::pair<int, int> >* alignment = 0,
               char gap = '-') const;

    /** Return if two sequences can be globally aligned.
    \note This method calls set_first() and set_second().
    */
    bool aligned(const std::string& first, const std::string& second);

private:
    mutable std::vector<int> matrix_;
    int gap_range_, max_errors_;
    const char* first_start_;
    const char* second_start_;
    int first_size_, second_size_;
    bool no_tail_;

    void cut_tail(int& first_last, int& second_last) const;

    void export_alignment(int row, int col,
                          std::string* first_str, std::string* second_str,
                          std::vector<std::pair<int, int> >* alignment,
                          char gap) const;

    int rows() const {
        return first_size_;
    }

    int cols() const {
        return second_size_;
    }

    int side() const;

    int row_size() const;

    int max_row() const;

    int min_col(int row) const;

    int max_col(int row) const;

    int& at(int row, int col) const;

    bool in(int row, int col) const;

    int substitution(int row, int col) const;

    void adjust_matrix_size();
};

}

#endif

